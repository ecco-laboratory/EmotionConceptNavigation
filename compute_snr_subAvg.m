[project_dir, fmri_data_dir, ~, theta_dir, ~, ~, ~, ~, spm_timing_dir, motion_regressors_dir, nuisance_regressors_dir, spm1stlevel_dir, ~, phi_dir, subjects] = set_up_dirs_constants();

smooth = 'unsmoothed'; 
data_dir = fullfile(project_dir, 'outputs', 'snr', smooth); 

files = dir(fullfile(data_dir, 'sub*.nii'));
% Create cell array of full file paths for all subject files
files = arrayfun(@(f) fullfile(f.folder, f.name), files, 'UniformOutput', false);

data = fmri_data(files); %voxel*subjects
%average across subjects (ignore nan)
avg_snr = nanmean(data.dat, 2); %1*voxel
%write to nifti
temp_nifti = fmri_data('/home/data/eccolab/MNS/data/fmri/nifti/derivatives/fmriprep-25.0.0/sub-0001/func/sub-0001_task-face_run-01_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii,1');
temp_nifti.dat = avg_snr;
temp_nifti.fullpath = fullfile(data_dir, 'subAvg_snr_avgRuns_wholeBrain.nii');
temp_nifti.write;