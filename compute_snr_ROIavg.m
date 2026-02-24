[project_dir, fmri_data_dir, ~, theta_dir, ~, ~, ~, ~, spm_timing_dir, motion_regressors_dir, nuisance_regressors_dir, spm1stlevel_dir, ~, phi_dir, subjects] = set_up_dirs_constants();

smooth = 'unsmoothed'; 
data_dir = fullfile(project_dir, 'outputs', 'snr', smooth); 

region_names = {'OFC2016ConstantinescuR5','vmPFC2019BaoR5','vmPFCcurrentStudyR5','vmPFC2016ConstantinescuR5','HC','ERC','alERC','pmERC'};
region_masks = {fullfile(project_dir, 'masks', 'OFC_2016Constantinescu_r5.nii'),...
                fullfile(project_dir, 'masks', 'vmPFC_2019Bao_r5.nii'),...
                fullfile(project_dir, 'masks', 'vmPFC_currentStudy_r5.nii'),...
                fullfile(project_dir, 'masks', 'vmPFC_2016Constantinescu_r5.nii'),...
                fullfile(project_dir, 'masks', 'HC_Julich.nii'),...
                fullfile(project_dir, 'masks', 'ERC_Julich.nii'),...
                fullfile(project_dir, 'masks', 'alEC_PRCpref_MNI.nii'),...
                fullfile(project_dir, 'masks', 'pmEC_PHCpref_MNI.nii')};

files = dir(fullfile(data_dir, 'sub*.nii'));

T = table();
row = 0;

for s = 1:numel(subjects)
    fname = fullfile(data_dir, sprintf('sub%s_snr_avgRuns_wholeBrain.nii', subjects{s}));
    data = fmri_data(fname);
    if ~isfile(fname), continue; end
    for r = 1:numel(region_names)
        region_name_fname = region_names{r};

        X = apply_mask(data, region_masks{r}).dat;
        vals = X(~isnan(X));

        row = row + 1;
        T.Subject(row,1)    = string(subjects{s});
        T.ROI(row,1)        = string(region_names{r});
        T.mean_tSNR(row,1)  = mean(vals);
        T.median_tSNR(row,1)= median(vals);
        T.nvox(row,1)       = numel(vals);
    end
end

%save to csv
%make a dir of data_dir, table
if ~exist(fullfile(data_dir, 'table'), 'dir')
    mkdir(fullfile(data_dir, 'table'))
end
if ~exist(fullfile(data_dir, 'table', 'snr_ROIavg.csv'), 'file')
    writetable(T, fullfile(data_dir, 'table', 'snr_ROIavg.csv'));
else
    %append to existing file
    writetable(T, fullfile(data_dir, 'table', 'snr_ROIavg.csv'), 'WriteMode', 'append');
end
