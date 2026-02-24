[project_dir, fmri_data_dir, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, subjects] = set_up_dirs_constants();

mean_signal_var_dir = fullfile(project_dir, 'outputs', 'mean_signal_var');
files = dir(fullfile(mean_signal_var_dir, 'all_subs', 'sub-*.nii'));
for f=1:length(files)
    P{f}=[files(f).folder filesep files(f).name];
end
dat=fmri_data(P);

x = mean(dat.dat')./std(dat.dat');
% template for writing nifti
temp_nifti = fmri_data('/home/data/eccolab/MNS/data/fmri/nifti/derivatives/fmriprep-25.0.0/sub-0001/func/sub-0001_task-face_run-01_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii,1');
temp_nifti.dat = x';
temp_nifti.fullpath = fullfile(mean_signal_var_dir, 'mean_signal_var_allRuns_across_subs.nii');
temp_nifti.write;
