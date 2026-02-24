[project_dir, fmri_data_dir, ~, theta_dir, ~, ~, ~, ~, spm_timing_dir, motion_regressors_dir, nuisance_regressors_dir, spm1stlevel_dir, ~, phi_dir, subjects] = set_up_dirs_constants();

TR = 1.25; % TR in seconds
discard_time = 7.5; % Time to discard in seconds
discard_volumes = discard_time / TR; % Number of volumes to discard 

modalities = {'word', 'face'};
runs = [1, 2]; % Two runs for each modality
smooth = 'unsmoothed'; 
output_dir = fullfile(project_dir, 'outputs', 'snr', smooth); if ~exist(output_dir, 'dir'), mkdir(output_dir); end

temp_nifti = fmri_data(fullfile(project_dir, 'data', 'fmri', 'nifti', 'derivatives', 'fmriprep-25.0.0', ['sub-', subjects{1}], 'func', ['sub-', subjects{1}, '_task-face_run-01_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii']));
temp_nifti.dat = temp_nifti.dat(:, 1);
for s = 1:length(subjects)
    subject = subjects{s};
    concat_runs = [];
    for m = 1:length(modalities)
        modality = modalities{m};
        for r = 1:length(runs)
            run = runs(r);
            disp(['Computing SNR for subject: ', subject, ' (', modality, ' run ', num2str(run), ')']);
            fname = fullfile(fmri_data_dir, ['sub-', subject], 'func', ['sub-', subject, '_task-', modality, '_run-0', num2str(run), '_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii']);
            

            data = fmri_data(fname).dat; %voxel*time 
            %discard from the first 'discard_volumes' volumes
            data = data(:, (discard_volumes+1):end);

            snr = mean(data') ./ std(data'); %1*voxel
            %delete to save memory
            clear data;
            concat_runs = [concat_runs; snr];
        end
    end
    %average across runs
    avg_snr = mean(concat_runs); %1*voxel
    snr_nifti = temp_nifti;
    snr_nifti.dat = avg_snr';
    snr_nifti.fullpath = fullfile(output_dir, ['sub', subject, '_snr_avgRuns_wholeBrain.nii']);
    %save nifti
    snr_nifti.write;
    
end