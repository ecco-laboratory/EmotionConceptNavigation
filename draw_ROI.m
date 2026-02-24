[project_dir, fmri_data_dir, ~, theta_dir, ~, ~, ~, ~, spm_timing_dir, motion_regressors_dir, nuisance_regressors_dir, spm1stlevel_dir, ~, phi_dir, subjects] = set_up_dirs_constants();
spm1stlevel_dir = [spm1stlevel_dir, '_noICA', '_wButton'];

addpath(fullfile(spm1stlevel_dir, ['sub', subjects{1}], 'singleTrial'))

clusters=sphere_mask('beta_0001.nii', [6, 44, -10], 5, fullfile(project_dir, 'masks', 'OFC_2016Constantinescu_r5.nii'));
V=spm_vol('/home/data/eccolab/MNS/outputs/spm_level1models_noICA_wButton/sub0008/singleTrial/beta_0001.nii');
V.fname = fullfile(project_dir, 'masks', 'OFC_2016Constantinescu_r5.nii');
m = clusters2mask(clusters, V.dim(1:3));
spm_write_vol(V, m);

clusters=sphere_mask('beta_0001.nii', [-8, 42, 0], 5, fullfile(project_dir, 'masks', 'vmPFC_2016Constantinescu_r5.nii'));
V=spm_vol('/home/data/eccolab/MNS/outputs/spm_level1models_noICA_wButton/sub0008/singleTrial/beta_0001.nii');
V.fname = fullfile(project_dir, 'masks', 'vmPFC_2016Constantinescu_r5.nii');
m = clusters2mask(clusters, V.dim(1:3));
spm_write_vol(V, m);


clusters=sphere_mask('beta_0001.nii', [6, 46, -10], 5, fullfile(project_dir, 'masks', 'vmPFC_2019Bao_r5.nii'));
V=spm_vol('/home/data/eccolab/MNS/outputs/spm_level1models_noICA_wButton/sub0008/singleTrial/beta_0001.nii');
V.fname = fullfile(project_dir, 'masks', 'vmPFC_2019Bao_r5.nii');
m = clusters2mask(clusters, V.dim(1:3));
spm_write_vol(V, m);

clusters=sphere_mask('beta_0001.nii', [2, 28, -20], 5, fullfile(project_dir, 'masks', 'vmPFC_2019Bao2_r5.nii'));
V=spm_vol('/home/data/eccolab/MNS/outputs/spm_level1models_noICA_wButton/sub0008/singleTrial/beta_0001.nii');
V.fname = fullfile(project_dir, 'masks', 'vmPFC_2019Bao2_r5.nii');
m = clusters2mask(clusters, V.dim(1:3));
spm_write_vol(V, m);