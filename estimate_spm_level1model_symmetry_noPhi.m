[project_dir, fmri_data_dir, ~, theta_dir, ~, ~, ~, ~, spm_timing_dir, motion_regressors_dir, nuisance_regressors_dir, spm1stlevel_dir, ~, ~, subjects] = set_up_dirs_constants();

spm1stlevel_dir = [spm1stlevel_dir, '_noICA', '_wButton'];
confounds_to_use = 'motion';%'motion_rejectedICA'; %'motion'
add_button = 'Button';%'';
theta_sources = {'Subavg', 'Subspec'};%,'Norm', 'Subspec'};%'Subavg'};%,'Subavg',  'Norm', 'Subspec'};
non_norm_theta_trials = 'includeNonNorm';%'excludeNonNorm';%'';
smooth = 'unsmoothed';%'unsmoothed';
if strcmp(smooth, 'smoothed'), smooth_data_prefix = 'smoothed_4mm_'; else smooth_data_prefix = ''; end
if strcmp(smooth, 'smoothed'), output_dir_base = spm1stlevel_dir; else output_dir_base = fullfile(spm1stlevel_dir, 'unsmoothed'); end
periodicity = {3};% {6, 4, 5, 7, 8};

% figure out which subjects this job should handle
%batch_size = 3;start_idx = (array_id - 1) * batch_size + 1; end_idx = min(start_idx + batch_size - 1, length(subjects)); subjects = subjects(start_idx:end_idx);
subjects = {subjects{array_id}};
disp(subjects)

modalities = {'word', 'face'};
runs = [1, 2]; % Two runs for each modality
TR = 1.25; % TR in seconds
discard_time = 7.5; % Time to discard in seconds
discard_volumes = discard_time / TR; % Number of volumes to discard 

for p = 1:length(periodicity)
    current_periodicity = periodicity{p};
    for t = 1:length(theta_sources)
        theta_source = theta_sources{t};
        for s = 1:length(subjects)
            subject = subjects{s};
            
            disp(['Processing subject: ', subject, ', periodicity: ', num2str(current_periodicity), ', theta source: ', theta_source]);
            
            if length(modalities) == 1
                subj_output_dir = fullfile(output_dir_base, ['sub', subject], 'symmetry_noPhi', non_norm_theta_trials, modalities{1}, ['periodicity', num2str(current_periodicity)], ['theta', theta_source]);
            else
                subj_output_dir = fullfile(output_dir_base, ['sub', subject], 'symmetry_noPhi', non_norm_theta_trials, ['periodicity', num2str(current_periodicity)], ['theta', theta_source]);
            end
            if ~exist(subj_output_dir, 'dir')
                mkdir(subj_output_dir);
            end
            
            matlabbatch = {};
            
            matlabbatch{1}.spm.stats.fmri_spec.dir = {subj_output_dir};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 36;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = round(36/2); 
            
            %% LOOP THROUGH MODALITIES AND RUNS TO ADD SESSIONS
            session_counter = 0;
            num_nuisance_regressors_all_sessions = [];
            num_conditions_all_sessions = [];
            face_sessions = [];
            word_sessions = [];
            cond_names_all_sessions = {};
            for m = 1:length(modalities)
                modality = modalities{m};
                for r = 1:length(runs)
                    run = runs(r);
                    
                    
                    %% Add image data for this session
                    func_files = spm_select('FPList', fullfile(fmri_data_dir, ['sub-', subject], 'func'), ...
                        ['^', smooth_data_prefix, 'sub-', subject, '_task-', modality, '_run-0', num2str(run), '_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.*\.nii$']);
                    if isempty(func_files)
                        warning(['No functional files found for subject ', subject, ', modality ', modality, ', run ', num2str(run)]);
                        continue;
                    else
                        session_counter = session_counter + 1;
                    end
                    V = spm_vol(func_files);
                    num_vols = length(V);
                    % skip the first 'discard_volumes' volumes 
                    func_files_cell = cell(num_vols - discard_volumes, 1);
                    for v = 1:(num_vols - discard_volumes)
                        actual_vol = v + discard_volumes;
                        func_files_cell{v} = [func_files, ',', num2str(actual_vol)];
                    end
                    
                    matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).scans = func_files_cell;

                    %% Add conditions for this session
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
                    %in timing file, there will be a variable 'names' that is a cell array with {'stimulus_name', 'rating_name', 'button_name'}
                    if run == 1
                        onsets = onsets_run1;      
                        durations = durations_run1; 
                    else
                        onsets = onsets_run2;
                        durations = durations_run2;
                    end
                    
                    % get theta for this subject and modality
                    theta_file = fullfile(theta_dir, ['sub', subject], ['sub', subject, '_', modality, '_thetas_run', num2str(run), '.mat']);
                    load(theta_file);
                    if strcmp(theta_source, 'Norm')
                        theta = thetas.Norm;
                    elseif strcmp(theta_source, 'Subspec')
                        theta = thetas.Subspec;
                        %if there's nan in theta, fill the nans with values from subavg
                        if any(isnan(theta))
                            theta(isnan(theta)) = thetas.Subavg(isnan(theta));
                        end
                    elseif strcmp(theta_source, 'Subavg')
                        theta = thetas.Subavg;
                        %get 'wrong' theta trials (if theta.Subspec is nan or its difference from theta.Subavg is greater than 10 degrees)
                        %notice the degree for both theta.Subspec and theta.Subavg is from 0 to 360 so 1 is different from 359 by 2 degree so can't just do abs(theta.Subspec - theta.Subavg) < 10
                        %so we need to consider the circular nature of the data
                        diff_theta = circular_diff(thetas.Subspec, thetas.Subavg);
                        %get wrong trial indecies
                        nonNorm_theta_trials = isnan(diff_theta) | diff_theta > 10;
                    end

                    cosTheta = cos(current_periodicity * deg2rad(theta));
                    sinTheta = sin(current_periodicity * deg2rad(theta));

                    if strcmp(non_norm_theta_trials, 'excludeNonNorm')
                        num_conditions = length(names) + 3; % 2 for cosTheta and sinTheta parametric modulators and 1 for wrong theta trials
                        num_conditions_all_sessions = [num_conditions_all_sessions, num_conditions];
                        cond_names_all_sessions{session_counter} = ['nonNormTheta', names(1), {'cosTheta', 'sinTheta'}, names(2:end)];
                        cond_idx = 1;
                        for cond_loc = 1:(length(names)+1)
                            orig_cond_name = names{cond_loc};
                            if strcmp(orig_cond_name, 'face') || strcmp(orig_cond_name, 'word')
                                %wrong theta trials
                                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).name = [orig_cond_name, '_nonNormTheta'];
                                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).onset = onsets{cond_idx}(nonNorm_theta_trials);
                                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).duration = durations{cond_idx}(nonNorm_theta_trials);
                                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).tmod = 0;
                                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).orth = 0;
                                %correct theta trials
                                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx+1).name = names{cond_loc};
                                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx+1).onset = onsets{cond_idx}(~nonNorm_theta_trials);
                                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx+1).duration = durations{cond_idx}(~nonNorm_theta_trials);
                                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx+1).tmod = 0;
                                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx+1).orth = 0;
                                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx+1).pmod = struct('name', 'cosTheta', 'param', cosTheta{~nonNorm_theta_trials}, 'poly', 1);
                                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx+2).pmod = struct('name', 'sinTheta', 'param', sinTheta{~nonNorm_theta_trials}, 'poly', 1);
                            else
                                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).name = names{cond_loc};
                                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).onset = onsets{cond_idx};
                                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).duration = durations{cond_idx};
                                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).tmod = 0;
                                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).orth = 0;
                                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).pmod = struct('name', {}, 'param', {}, 'poly', {});
                                cond_idx = cond_idx + 1;
                            end
                        end

                    else 
                        num_conditions = length(names) + 2; % 2 for cosTheta and sinTheta parametric modulators
                        num_conditions_all_sessions = [num_conditions_all_sessions, num_conditions];
                        %add cos and sin theta to the names (should be after the 1st condition in the names cell array and before the 2nd condition)
                        cond_names_all_sessions{session_counter} = [names(1), {'cosTheta', 'sinTheta'}, names(2:end)];

                        for cond_idx = 1:length(names)
                            matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).name = names{cond_idx};
                            matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).onset = onsets{cond_idx};
                            matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).duration = durations{cond_idx};
                            matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).tmod = 0;
                            matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).orth = 0;
                            
                            % if this is a stimulus condition (face or word)
                            orig_cond_name = names{cond_idx};
                            if strcmp(orig_cond_name, 'face') || strcmp(orig_cond_name, 'word')
                                % Get the cosine and sine of the theta
                                
                                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).pmod(1) = struct('name', 'cosTheta', 'param', cosTheta, 'poly', 1);
                                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).pmod(2) = struct('name', 'sinTheta', 'param', sinTheta, 'poly', 1);
                            else
                                matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).cond(cond_idx).pmod = struct('name', {}, 'param', {}, 'poly', {});
                            end
                        end
                    end

                    %% Add nuisance regressors for this session
                    if strcmp(confounds_to_use, 'motion_rejectedICA')
                        confounds_file = fullfile(nuisance_regressors_dir, ['sub', subject], ['sub', subject, '_motion_tedanaRejectedICA_confounds_', modality, '_run', num2str(run), '.txt']);
                    else
                        confounds_file = fullfile(motion_regressors_dir, ['sub', subject], ['sub', subject, '_motionRegressors_', modality, '_run', num2str(run), '.txt']);
                    end
                    if exist(confounds_file, 'file')
                        matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).multi = {''};
                        matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).regress = struct('name', {}, 'val', {});
                        matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).multi_reg = {confounds_file};
                        num_nuisance_regressors = size(readtable(confounds_file), 2);
                        num_nuisance_regressors_all_sessions = [num_nuisance_regressors_all_sessions, num_nuisance_regressors];
                    else
                        error(['Confounds file not found: ', confounds_file]);
                    end
                    
                    matlabbatch{1}.spm.stats.fmri_spec.sess(session_counter).hpf = 128;
                    
                    if strcmp(modality, 'face')
                        face_sessions = [face_sessions, session_counter];
                    else
                        word_sessions = [word_sessions, session_counter];
                    end
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
            
            %% Model estimation
            matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(subj_output_dir, 'SPM.mat')};
            matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
            
            %% Contrast specification 
            matlabbatch{3}.spm.stats.con.spmmat = {fullfile(subj_output_dir, 'SPM.mat')};


            total_num_regressors = sum(num_conditions_all_sessions) + sum(num_nuisance_regressors_all_sessions);
            total_num_regressors_per_session = num_conditions_all_sessions+num_nuisance_regressors_all_sessions;%array of num of regressors per session
            contrast_count = 1;
            
            total_sessions = length(face_sessions) + length(word_sessions);
            regressor_map = cell(1, total_sessions);    
            for sess = 1:total_sessions
                base_idx = sum(total_num_regressors_per_session(1:sess-1));
                cond_names = cond_names_all_sessions{sess};
                session_struct = struct();
                for cond_idx = 1:length(cond_names)
                    cond_name = matlab.lang.makeValidName(cond_names{cond_idx}); % ensure valid field name
                    session_struct.(cond_name) = base_idx + cond_idx;
                end
                regressor_map{sess} = session_struct;
            end

            % 1. Contrast: cosTheta only
            cosTheta_contrast = zeros(1, total_num_regressors);
            if ~isempty(face_sessions)
                for sess = face_sessions
                    cosTheta_contrast(regressor_map{sess}.cosTheta) = 1/(length(face_sessions)+length(word_sessions));
                end
            end
            if ~isempty(word_sessions)
                for sess = word_sessions
                    cosTheta_contrast(regressor_map{sess}.cosTheta) = 1/(length(face_sessions)+length(word_sessions));
                end
            end
            matlabbatch{3}.spm.stats.con.consess{contrast_count}.tcon.name = 'cosTheta';
            matlabbatch{3}.spm.stats.con.consess{contrast_count}.tcon.weights = cosTheta_contrast;
            matlabbatch{3}.spm.stats.con.consess{contrast_count}.tcon.sessrep = 'none';
            contrast_count = contrast_count + 1;

            % 2. Contrast: sinTheta only
            sinTheta_contrast = zeros(1, total_num_regressors);
            if ~isempty(face_sessions)
                for sess = face_sessions
                    sinTheta_contrast(regressor_map{sess}.sinTheta) = 1/(length(face_sessions)+length(word_sessions));
                end
            end
            if ~isempty(word_sessions)
                for sess = word_sessions
                    sinTheta_contrast(regressor_map{sess}.sinTheta) = 1/(length(face_sessions)+length(word_sessions));
                end
            end
            matlabbatch{3}.spm.stats.con.consess{contrast_count}.tcon.name = 'sinTheta';
            matlabbatch{3}.spm.stats.con.consess{contrast_count}.tcon.weights = sinTheta_contrast;
            matlabbatch{3}.spm.stats.con.consess{contrast_count}.tcon.sessrep = 'none';
            contrast_count = contrast_count + 1;
            
            % 3. Contrast: ftest (two rows, one cosTheta, one sinTheta)
            ftest_contrast = [cosTheta_contrast; sinTheta_contrast];
            matlabbatch{3}.spm.stats.con.consess{contrast_count}.fcon.name = 'ftest_CosSinTheta';
            matlabbatch{3}.spm.stats.con.consess{contrast_count}.fcon.weights = ftest_contrast;
            matlabbatch{3}.spm.stats.con.consess{contrast_count}.fcon.sessrep = 'none';
            contrast_count = contrast_count + 1;

            % 4. Contrast: face sinTheta
            if ~isempty(face_sessions) && ~isempty(word_sessions)
                face_sinTheta_contrast = zeros(1, total_num_regressors);
                face_cosTheta_contrast = zeros(1, total_num_regressors);
                for sess = face_sessions
                    face_sinTheta_contrast(regressor_map{sess}.sinTheta) = 1/(length(face_sessions));
                    face_cosTheta_contrast(regressor_map{sess}.cosTheta) = 1/(length(face_sessions));
                end
                matlabbatch{3}.spm.stats.con.consess{contrast_count}.tcon.name = 'face_sinTheta';
                matlabbatch{3}.spm.stats.con.consess{contrast_count}.tcon.weights = face_sinTheta_contrast;
                matlabbatch{3}.spm.stats.con.consess{contrast_count}.tcon.sessrep = 'none';
                contrast_count = contrast_count + 1;

            % 5. Contrast: face cosTheta
                matlabbatch{3}.spm.stats.con.consess{contrast_count}.tcon.name = 'face_cosTheta';
                matlabbatch{3}.spm.stats.con.consess{contrast_count}.tcon.weights = face_cosTheta_contrast;
                matlabbatch{3}.spm.stats.con.consess{contrast_count}.tcon.sessrep = 'none';
                contrast_count = contrast_count + 1;

            % 6. Contrast: face ftest
                face_ftest_contrast = [face_sinTheta_contrast; face_cosTheta_contrast];
                matlabbatch{3}.spm.stats.con.consess{contrast_count}.fcon.name = 'face_ftest_SinCosTheta';
                matlabbatch{3}.spm.stats.con.consess{contrast_count}.fcon.weights = face_ftest_contrast;
                matlabbatch{3}.spm.stats.con.consess{contrast_count}.fcon.sessrep = 'none';
                contrast_count = contrast_count + 1;
            
            % 7. Contrast: word sinTheta
                word_sinTheta_contrast = zeros(1, total_num_regressors);
                word_cosTheta_contrast = zeros(1, total_num_regressors);
                for sess = word_sessions
                    word_sinTheta_contrast(regressor_map{sess}.sinTheta) = 1/(length(word_sessions));
                    word_cosTheta_contrast(regressor_map{sess}.cosTheta) = 1/(length(word_sessions));
                end
                matlabbatch{3}.spm.stats.con.consess{contrast_count}.tcon.name = 'word_sinTheta';
                matlabbatch{3}.spm.stats.con.consess{contrast_count}.tcon.weights = word_sinTheta_contrast;
                matlabbatch{3}.spm.stats.con.consess{contrast_count}.tcon.sessrep = 'none';
                contrast_count = contrast_count + 1;

            % 8. Contrast: word cosTheta
                matlabbatch{3}.spm.stats.con.consess{contrast_count}.tcon.name = 'word_cosTheta';
                matlabbatch{3}.spm.stats.con.consess{contrast_count}.tcon.weights = word_cosTheta_contrast;
                matlabbatch{3}.spm.stats.con.consess{contrast_count}.tcon.sessrep = 'none';
                contrast_count = contrast_count + 1;

            % 9. Contrast: word ftest
                word_ftest_contrast = [word_sinTheta_contrast; word_cosTheta_contrast];
                matlabbatch{3}.spm.stats.con.consess{contrast_count}.fcon.name = 'word_ftest_SinCosTheta';
                matlabbatch{3}.spm.stats.con.consess{contrast_count}.fcon.weights = word_ftest_contrast;
                matlabbatch{3}.spm.stats.con.consess{contrast_count}.fcon.sessrep = 'none';
                contrast_count = contrast_count + 1;
            end
            
            matlabbatch{3}.spm.stats.con.delete = 0;
            
            try
                spm_jobman('run', matlabbatch);
                disp(['Successfully completed processing for subject ', subject]);
                
            catch ME
                warning(['Error processing subject ', subject, ': ', ME.message]);
            end
        end
    end
end
