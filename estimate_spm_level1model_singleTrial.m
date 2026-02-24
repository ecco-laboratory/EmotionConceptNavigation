clear; close all; clc;
[project_dir, fmri_data_dir, ~, ~, ~, ~, ~, ~, spm_timing_dir, motion_regressors_dir, nuisance_regressors_dir, spm1stlevel_dir, ~, ~, subjects] = set_up_dirs_constants();

spm1stlevel_dir = [spm1stlevel_dir, '_noICA', '_wButton']; %[spm1stlevel_dir, '_noICA', '_wButton'];
confounds_to_use = 'motion'; %'motion'; %'motion_rejectedICA'; %'motion'
add_button = 'Button';%'';
modalities = {'word', 'face'};
runs = [1, 2]; % Two runs for each modality
TR = 1.25; % TR in seconds
discard_time = 7.5; % Time to discard in seconds
discard_volumes = discard_time / TR; % Number of volumes to discard 
smooth = 'unsmoothed';%'smoothed';%'';
if strcmp(smooth, 'smoothed'), smooth_data_prefix = 'smoothed_4mm_'; else smooth_data_prefix = ''; end

spm('defaults','fmri');
spm_jobman('initcfg');

for s = 1:length(subjects)
    subject = subjects{s};
    subj_output_dir = fullfile(spm1stlevel_dir, smooth, ['sub', subject], 'singleTrial');
    if ~exist(subj_output_dir, 'dir')
        mkdir(subj_output_dir);
    end
    
    matlabbatch = {};
    matlabbatch{1}.spm.stats.fmri_spec.dir = {subj_output_dir};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 36;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = round(36/2);
    

    session_counter = 0;
    %% LOOP THROUGH MODALITIES AND RUNS TO ADD SESSIONS
    for m = 1:length(modalities)
        modality = modalities{m};
        
        for r = 1:length(runs)
            run = runs(r);
            session_counter = session_counter + 1;
            
            disp(['Processing subject: ', subject, ', modality: ', modality, ', run: ', num2str(run), ' (session ', num2str(session_counter), ')']);
            
            func_files = spm_select('FPList', fullfile(fmri_data_dir, ['sub-', subject], 'func'), ...
                ['^', smooth_data_prefix, 'sub-', subject, '_task-', modality, '_run-0', num2str(run), '_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.*\.nii$']);
            if isempty(func_files)
                warning(['No image data found for subject ', subject, ', modality ', modality, ', run ', num2str(run)]);
                continue;
            end
            V = spm_vol(func_files); %number of volumes
            num_vols = length(V);
            %skip the first 'discard_volumes' volumes 
            func_files_cell = cell(num_vols - discard_volumes, 1);
            for v = 1:(num_vols - discard_volumes)
                actual_vol = v + discard_volumes;
                func_files_cell{v} = [func_files, ',', num2str(actual_vol)];
            end
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).scans = func_files_cell;
            
            % Load timing files 
            if strcmp(modality, 'word')
                timing_file = fullfile(spm_timing_dir, ['sub', subject], ['sub', subject, '_', modality, 'RatingContext', add_button, '_run', num2str(run), '.mat']);
            else
                timing_file = fullfile(spm_timing_dir, ['sub', subject], ['sub', subject, '_', modality, 'Rating', add_button, '_run', num2str(run), '.mat']);
            end
            if ~exist(timing_file, 'file')
                warning(['Timing file not found: ', timing_file]);
                continue;
            end
            load(timing_file);
            %in timing file, there is also be a variable 'names' that is a cell array with {'stimulus_name', 'rating_name', 'button_name'}
            if run == 1
                onsets = onsets_run1;      % cell array with {stim_onsets, rating_onsets, button_onsets}
                durations = durations_run1; % cell array with {stim_durations, rating_durations, button_durations}
            else
                onsets = onsets_run2;
                durations = durations_run2;
            end
            stim_onsets = onsets{1};    
            stim_durations = durations{1}; 

            %stimulus single trial conditions 
            cond_idx = 0;
            for trial = 1:length(stim_onsets)
                cond_idx = cond_idx + 1;
                
                trial_name = sprintf('%s_run%d_trial%03d', modality, run, trial);
                % Add this single trial as a separate condition
                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).name = trial_name;
                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).onset = stim_onsets(trial);
                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).duration = stim_durations(trial);
                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).tmod = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).pmod = struct('name', {}, 'param', {}, 'poly', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).orth = 1;
            end
            
            % Add rating and button press conditions (with modality prefixes to distinguish them)
            for other_cond_idx = 2:length(names)
                cond_idx = cond_idx + 1;
                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).name = names{other_cond_idx};
                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).onset = onsets{other_cond_idx};
                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).duration = durations{other_cond_idx};
                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).tmod = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).pmod = struct('name', {}, 'param', {}, 'poly', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).orth = 1;
            end
            
            % Load motion regressors 
            if strcmp(confounds_to_use, 'motion_rejectedICA')
                confounds_file = fullfile(nuisance_regressors_dir, ['sub', subject], ['sub', subject, '_motion_tedanaRejectedICA_confounds_', modality, '_run', num2str(run), '.txt']);
            else
                confounds_file = fullfile(motion_regressors_dir, ['sub', subject], ['sub', subject, '_motionRegressors_', modality, '_run', num2str(run), '.txt']);
            end
            if exist(confounds_file, 'file')
                % Add 6 motion regressors
                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).multi = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).regress = struct('name', {}, 'val', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).multi_reg = {confounds_file};
            else
                error(['Confounds file not found: ', confounds_file]);
            end
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).hpf = 128;

        end
    end
    
    %% MODEL ESTIMATION SETTINGS %%%%%
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = -Inf;

    spm_path = spm('Dir');
    icbm_mask = fullfile(spm_path, 'tpm', 'mask_ICV.nii');
    matlabbatch{1}.spm.stats.fmri_spec.mask = {icbm_mask};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(subj_output_dir, 'SPM.mat')};
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
    try
        spm_jobman('run', matlabbatch);  
        disp(['Successfully completed single trial analysis for subject ', subject]);
    catch ME
        warning(['Error processing subject ', subject, ': ', ME.message]);
    end
end