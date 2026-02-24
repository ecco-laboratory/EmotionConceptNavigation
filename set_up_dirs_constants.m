function [project_dir, fmri_data_dir, psychopy_csv_dir, theta_dir, category_dir, vaDiff_dir, avgVA_dir, distance_dir, spm_timing_dir, motion_regressors_dir, nuisance_regressors_dir, spm1stlevel_dir, spm2ndlevel_dir, phi_dir, subjects] = set_up_dirs_constants()
    addpath('/home/data/eccolab/Code/GitHub/spm12');
    addpath(genpath('/home/data/eccolab/Code/GitHub/CanlabCore'))
    addpath(genpath('/home/data/eccolab/Code/GitHub/Neuroimaging_Pattern_Masks/'))
    addpath(genpath('/home/data/eccolab/MNS/toolbox'))

    project_dir = fileparts(mfilename('fullpath'));
    fmri_data_dir = fullfile(project_dir, 'data', 'fmri', 'nifti', 'derivatives', 'fmriprep-25.0.0');
    psychopy_csv_dir = fullfile(project_dir, 'data', 'beh');
    theta_dir = fullfile(project_dir, 'data', 'beh', 'theta');
    category_dir = fullfile(project_dir, 'data', 'beh', 'category');
    vaDiff_dir = fullfile(project_dir, 'data', 'beh', 'vaDiff');
    avgVA_dir = fullfile(project_dir, 'data', 'beh', 'avgVA');
    distance_dir = fullfile(project_dir, 'data', 'beh', 'distance');
    spm_timing_dir = fullfile(project_dir, 'outputs', 'spm_condition_timing');
    motion_regressors_dir = fullfile(project_dir, 'outputs', 'motion_regressors');
    nuisance_regressors_dir = fullfile(project_dir, 'outputs', 'nuisance_regressors');
    spm1stlevel_dir = fullfile(project_dir, 'outputs', 'spm_level1models');
    spm2ndlevel_dir = fullfile(project_dir, 'outputs', 'spm_level2models');
    phi_dir = fullfile(project_dir, 'outputs', 'phi');
    subjects = {'0019','0017','0011', '0021','0008','0014', '0015', '0016','0018', '0013','0012', '0010', '0009', '0007','0004', '0003', '0002', '0001'};

end


