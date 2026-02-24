[project_dir, fmri_data_dir, ~, theta_dir, ~, ~, ~, ~, spm_timing_dir, motion_regressors_dir, nuisance_regressors_dir, spm1stlevel_dir, ~, phi_dir, subjects] = set_up_dirs_constants();

data_dir = fullfile(project_dir, 'outputs', 'singleTrialBetaAnalysis', 'noICA', 'incl_all_subs_trials', 'onoffGridcontrast_multivariate_searchlight', 'phi_current', 'Subavg', 'periodicity6');
output_dir = fullfile(data_dir, 'splitHalfCV_avg', 'phi_voxelComponentAverage','singleTrialBeta','nifti');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
cross_validations = {'xModalityRun', 'xModality'};
region_names = {'wholebrain'};%{'HC', 'ERC', 'vmPFCGlasser'};



%average across folds for each subject
for s = 1:length(subjects)
    current_subject = subjects{s};
    for r = 1:length(region_names)
        current_region_name = region_names{r};
        %find all files under data_dir for current_subject ('contrast_value_sub<current_subject>_*_<region>.nii')
        files = dir(fullfile(data_dir, '**', sprintf('contrast_value_sub%s_*_%s.nii', current_subject, current_region_name)));
        paths = fullfile({files.folder}, {files.name});

        data = fmri_data(paths); %voxels * folds
        %average across folds
        data_avg = mean(data.dat, 2); %voxels * 1
        %save to nifti file
        temp_nifti = fmri_data(paths{1});
        temp_nifti.dat = data_avg;
        temp_nifti.fullpath = fullfile(output_dir, sprintf('contrast_value_sub%s_%s.nii', current_subject, current_region_name));
        temp_nifti.write;
    end
end

%group level analysis
for r = 1:length(region_names)
    current_region_name = region_names{r};
    %find all files under output_dir for current_region_name ('contrast_value_sub*_*.nii')
    files = dir(fullfile(output_dir, sprintf('contrast_value_sub*_%s.nii', current_region_name)));
    paths = fullfile({files.folder}, {files.name});
    data = fmri_data(paths); %voxels * subjects


    %now get cohen's d across subjects
    data_avg = mean(data.dat, 2); %voxels * 1
    data_std = std(data.dat, 0, 2); %voxels * 1
    cohen_d = data_avg ./ data_std;
    %save to nifti file
    temp_nifti = fmri_data(paths{1});
    temp_nifti.dat = cohen_d;
    temp_nifti.fullpath = fullfile(output_dir, sprintf('cohen_d_contrastValue_%s.nii', current_region_name));
    temp_nifti.write;

    %get t-statistic


    %{
    data_thre = threshold(ttest(data), .05, 'FDR');
    data_thre.fullpath = fullfile(output_dir, sprintf('tstat_contrastValue_FDR05_%s.nii', current_region_name));
    data_thre.write;
    data_thre = threshold(ttest(data), .05, 'unc');
    data_thre.fullpath = fullfile(output_dir, sprintf('tstat_contrastValue_UNC05_%s.nii', current_region_name));
    data_thre.write;
    %}
end

% manual t-test
for r = 1:length(region_names)

    current_region_name = region_names{r};
    files = dir(fullfile(output_dir, sprintf('contrast_value_sub*_%s.nii', current_region_name)));

    paths = fullfile({files.folder}, {files.name});
    data  = fmri_data(paths);

    X = data.dat; N = size(X, 2); df = N - 1;

    mu = mean(X, 2); sd = std(X, 0, 2); se = sd ./ sqrt(N); tvals = mu ./ se;
    
    % One-tailed p-values 
    pvals = 1 - tcdf(tvals, df);

    % Uncorrected threshold (p < .05, one-tailed)
    t_unc = tvals; t_unc(pvals > 0.05) = 0; 
    t_unc01 = tvals; t_unc01(pvals > 0.01) = 0;

    % FDR
    [p_sorted, sort_idx] = sort(pvals);
    V = length(pvals);
    q = 0.05;

    thresh_line = (1:V)'/V * q;
    below = p_sorted <= thresh_line;

    if any(below)
        max_idx = find(below, 1, 'last');
        p_thresh = p_sorted(max_idx);
        fdr_mask = pvals <= p_thresh;
    else
        fdr_mask = false(size(pvals));
    end

    t_fdr = tvals;
    t_fdr(~fdr_mask) = 0;

    % save maps
    tmp = fmri_data(paths{1});

    % Raw t
    tmp.dat = tvals;
    tmp.fullpath = fullfile(output_dir, sprintf('tstat_manual_pos_%s.nii', current_region_name));
    tmp.write;

    % Uncorrected
    tmp.dat = t_unc;
    tmp.fullpath = fullfile(output_dir, sprintf('tstat_manual_pos_UNC05_%s.nii', current_region_name));
    tmp.write;
    tmp.dat = t_unc01;
    tmp.fullpath = fullfile(output_dir, sprintf('tstat_manual_pos_UNC01_%s.nii', current_region_name));
    tmp.write;

    % FDR
    tmp.dat = t_fdr;
    tmp.fullpath = fullfile(output_dir, sprintf('tstat_manual_pos_FDR05_%s.nii', current_region_name));
    tmp.write;

    fprintf('Finished region: %s (N=%d)\n', current_region_name, N);

end

d = fmri_data(fullfile(output_dir, 'tstat_manual_pos_UNC01_wholebrain.nii'));
k_thresh = 20;
d2 = threshold(d, [0, Inf], 'raw-between', 'k', k_thresh);
d2.fullpath = fullfile(output_dir, sprintf('tstat_manual_pos_UNC01_k%d_wholebrain.nii', k_thresh));
d2.write;
%get peak coordinate table
d = fmri_data(fullfile(output_dir, sprintf('tstat_manual_pos_UNC01_k%d_wholebrain.nii', k_thresh)));
cl = region(d);

ncl = length(cl);
peak_t = nan(ncl,1);
peak_xyzmm = nan(ncl,3);
peak_xyz = nan(ncl,3);
k = nan(ncl,1);
M = nan(ncl,1);

for i = 1:ncl
    peak_t(i) = max(cl(i).Z);
    peak_xyzmm(i,:) = cl(i).mm_center;
    peak_xyz(i,:) = cl(i).center;
    k(i) = cl(i).numVox;
end

T = table( ...
    k, ...
    peak_t, ...
    peak_xyz(:,1), ...
    peak_xyz(:,2), ...
    peak_xyz(:,3), ...
    peak_xyzmm(:,1), ...
    peak_xyzmm(:,2), ...
    peak_xyzmm(:,3), ...
    'VariableNames', {'k_vox','peak_t','x','y','z','xmm','ymm','zmm'});


writetable(T, fullfile(output_dir, sprintf('peak_coords_UNC01_k%d_wholebrain.csv', k_thresh)));