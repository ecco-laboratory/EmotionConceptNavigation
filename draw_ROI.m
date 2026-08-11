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

clusters=sphere_mask('beta_0001.nii', [0, -32, 28], 5, fullfile(project_dir, 'masks', 'PCC_2016Constantinescu_r5.nii'));
V=spm_vol('/home/data/eccolab/MNS/outputs/spm_level1models_noICA_wButton/sub0008/singleTrial/beta_0001.nii');
V.fname = fullfile(project_dir, 'masks', 'PCC_2016Constantinescu_r5.nii');
m = clusters2mask(clusters, V.dim(1:3));
spm_write_vol(V, m);

clusters=sphere_mask('beta_0001.nii', [6, -52, 24], 5, fullfile(project_dir, 'masks', 'RSC_2016Constantinescu_r5.nii'));
V=spm_vol('/home/data/eccolab/MNS/outputs/spm_level1models_noICA_wButton/sub0008/singleTrial/beta_0001.nii');
V.fname = fullfile(project_dir, 'masks', 'RSC_2016Constantinescu_r5.nii');
m = clusters2mask(clusters, V.dim(1:3));
spm_write_vol(V, m);

clusters=sphere_mask('beta_0001.nii', [30, -62, 28], 5, fullfile(project_dir, 'masks', 'LPC_2016Constantinescu_r5.nii'));
V=spm_vol('/home/data/eccolab/MNS/outputs/spm_level1models_noICA_wButton/sub0008/singleTrial/beta_0001.nii');
V.fname = fullfile(project_dir, 'masks', 'LPC_2016Constantinescu_r5.nii');
m = clusters2mask(clusters, V.dim(1:3));
spm_write_vol(V, m);

clusters=sphere_mask('beta_0001.nii', [52, -42, 40], 5, fullfile(project_dir, 'masks', 'TPJ_2016Constantinescu_r5.nii'));
V=spm_vol('/home/data/eccolab/MNS/outputs/spm_level1models_noICA_wButton/sub0008/singleTrial/beta_0001.nii');
V.fname = fullfile(project_dir, 'masks', 'TPJ_2016Constantinescu_r5.nii');
m = clusters2mask(clusters, V.dim(1:3));
spm_write_vol(V, m);

for r = [11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    clusters=sphere_mask('beta_0001.nii', [6, 44, -10], r, fullfile(project_dir, 'masks', ['OFC_2016Constantinescu_r', num2str(r), '.nii']));
    V=spm_vol('/home/data/eccolab/MNS/outputs/spm_level1models_noICA_wButton/sub0008/singleTrial/beta_0001.nii');
    V.fname = fullfile(project_dir, 'masks', ['OFC_2016Constantinescu_r', num2str(r), '.nii']);
    m = clusters2mask(clusters, V.dim(1:3));
    spm_write_vol(V, m);
end

clusters=sphere_mask('beta_0001.nii', [-3, 63, 15], 5, fullfile(project_dir, 'masks', 'mPFC_2010Doeller_r5.nii'));
V=spm_vol('/home/data/eccolab/MNS/outputs/spm_level1models_noICA_wButton/sub0008/singleTrial/beta_0001.nii');
V.fname = fullfile(project_dir, 'masks', 'mPFC_2010Doeller_r5.nii');
m = clusters2mask(clusters, V.dim(1:3));
spm_write_vol(V, m);

clusters=sphere_mask('beta_0001.nii', [-54, 9, -30], 5, fullfile(project_dir, 'masks', 'LTCleft_2010Doeller_r5.nii'));
V=spm_vol('/home/data/eccolab/MNS/outputs/spm_level1models_noICA_wButton/sub0008/singleTrial/beta_0001.nii');
V.fname = fullfile(project_dir, 'masks', 'LTCleft_2010Doeller_r5.nii');
m = clusters2mask(clusters, V.dim(1:3));
spm_write_vol(V, m);

clusters=sphere_mask('beta_0001.nii', [42, 15, -36], 5, fullfile(project_dir, 'masks', 'LTCright_2010Doeller_r5.nii'));
V=spm_vol('/home/data/eccolab/MNS/outputs/spm_level1models_noICA_wButton/sub0008/singleTrial/beta_0001.nii');
V.fname = fullfile(project_dir, 'masks', 'LTCright_2010Doeller_r5.nii');
m = clusters2mask(clusters, V.dim(1:3));
spm_write_vol(V, m);

clusters=sphere_mask('beta_0001.nii', [-18, -54, 45], 5, fullfile(project_dir, 'masks', 'PPC_2010Doeller_r5.nii'));
V=spm_vol('/home/data/eccolab/MNS/outputs/spm_level1models_noICA_wButton/sub0008/singleTrial/beta_0001.nii');
V.fname = fullfile(project_dir, 'masks', 'PPC_2010Doeller_r5.nii');
m = clusters2mask(clusters, V.dim(1:3));
spm_write_vol(V, m);

clusters=sphere_mask('beta_0001.nii', [-0.5, 51.5, -16.5], 5, fullfile(project_dir, 'masks', 'vmPFC_currUnivariate_r5.nii'));
V=spm_vol('/home/data/eccolab/MNS/outputs/spm_level1models_noICA_wButton/sub0008/singleTrial/beta_0001.nii');
V.fname = fullfile(project_dir, 'masks', 'vmPFC_currUnivariate_r5.nii');
m = clusters2mask(clusters, V.dim(1:3));
spm_write_vol(V, m);


coords = [
    -0.5, 63.5,  -6.5
    -0.5, 49.5, -10.5
    -0.5, 53.5,  -6.5
     3.5, 61.5, -18.5
     3.5, 55.5,   3.5
     3.5, 57.5,  -4.5
   -10.5, 59.5,  -0.5
];


for i = [2, 3, 4, 5, 6, 7, 8]
    mask_name = sprintf('vmPFC_currUnivariate%d_r5.nii', i);
    clusters=sphere_mask('beta_0001.nii', coords(i-1, :), 5, fullfile(project_dir, 'masks', mask_name));
    V=spm_vol('/home/data/eccolab/MNS/outputs/spm_level1models_noICA_wButton/sub0008/singleTrial/beta_0001.nii');
    V.fname = fullfile(project_dir, 'masks', mask_name);
    m = clusters2mask(clusters, V.dim(1:3));
    spm_write_vol(V, m);
end