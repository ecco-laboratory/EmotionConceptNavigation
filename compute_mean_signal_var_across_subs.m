[project_dir, fmri_data_dir, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, subjects] = set_up_dirs_constants();

TR = 1.25;
discard_time = 7.5;
discard_volumes = round(discard_time / TR);

modalities = {'word', 'face'};
runs = [1, 2];

output_dir = fullfile(project_dir, 'outputs', 'mean_signal_var', 'all_subs');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% template for writing nifti
temp_nifti = fmri_data('/home/data/eccolab/MNS/data/fmri/nifti/derivatives/fmriprep-25.0.0/sub-0001/func/sub-0001_task-face_run-01_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii');
temp_nifti.dat = temp_nifti.dat(:,1);

for s = 1:length(subjects)
    subject = subjects{s};

    all_data = [];

    for m = 1:length(modalities)
        modality = modalities{m};
        for r = 1:length(runs)
            run = runs(r);

            fname = fullfile(fmri_data_dir, ['sub-', subject], 'func', ...
                ['sub-', subject, '_task-', modality, '_run-0', num2str(run), ...
                 '_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii']);

            disp(['Loading ', fname])
            data = fmri_data(fname).dat;   % voxel × time
            data = data(:, discard_volumes+1:end);

            all_data = [all_data, data];
            clear data
        end
    end

    % voxelwise stats across all runs
    curr_temp_nifti = temp_nifti;
    curr_temp_nifti.dat = nanmean(all_data, 2);
    curr_temp_nifti.fullpath = fullfile(output_dir, ['sub-', subject, '_mean_signal_allRuns.nii']);
    curr_temp_nifti.write;
end
