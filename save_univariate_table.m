[project_dir, fmri_data_dir, ~, theta_dir, ~, ~, ~, ~, spm_timing_dir, motion_regressors_dir, nuisance_regressors_dir, spm1stlevel_dir, ~, phi_dir, subjects] = set_up_dirs_constants();

t_map_dir = fullfile(project_dir, 'outputs', 'spm_level2models_noICA_wButton', 'symmetry_noPhi_Z_spm', 'includeNonNorm', 'SPM', 'periodicity6', 'thetaSubavg', 'ftest_CosSinTheta');
tvals = fmri_data(fullfile(t_map_dir, 'spmT_0001.nii')).dat;
df = 18-1;
pvals = 1 - tcdf(tvals, df);

 % Uncorrected threshold (p < .05, one-tailed)
 %t_unc = tvals; t_unc(pvals > 0.05) = 0; 
 %t_unc01 = tvals; t_unc01(pvals > 0.01) = 0;
 t_unc0001 = tvals; t_unc0001(pvals > 0.001) = 0; 

 % FDR
 [p_sorted, sort_idx] = sort(pvals);
 V = length(pvals);
 q = 0.01;

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
 tmp = fmri_data(fullfile(t_map_dir, 'spmT_0001.nii'));

  % FDR
 tmp.dat = t_fdr;
 k_thresh = 10;
 tmp = threshold(tmp, [0, Inf], 'raw-between', 'k', k_thresh);
 tmp.fullpath = fullfile(t_map_dir, 'thres_results', sprintf('tstat_manual_pos_FDR01_k%d_wholebrain.nii', k_thresh));
 tmp.write;


 % Uncorrected
 %tmp.dat = t_unc;
 %tmp.fullpath = fullfile(t_map_dir, 'thres_results', 'tstat_manual_pos_UNC05_wholebrain.nii');
 %tmp.write;
 tmp.dat = t_unc0001;
 tmp.fullpath = fullfile(t_map_dir, 'thres_results', 'tstat_manual_pos_UNC0001_wholebrain.nii');
 tmp.write;




d = fmri_data(fullfile(t_map_dir, 'thres_results', 'tstat_manual_pos_UNC0001_wholebrain.nii'));
k_thresh = 0;
d2 = threshold(d, [0, Inf], 'raw-between', 'k', k_thresh);
d2.fullpath = fullfile(t_map_dir, 'thres_results', sprintf('tstat_manual_pos_UNC0001_k%d_wholebrain.nii', k_thresh));
d2.write;
%get peak coordinate table
d = fmri_data(fullfile(t_map_dir, 'thres_results', sprintf('tstat_manual_pos_UNC0001_k%d_wholebrain.nii', k_thresh)));
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


writetable(T, fullfile(t_map_dir, 'thres_results', sprintf('peak_coords_UNC0001_k%d_wholebrain.csv', k_thresh)));



x=fmri_data('/home/data/eccolab/MNS/outputs/spm_level2models_noICA_wButton/symmetry_noPhi_Z_spm/includeNonNorm/SPM/periodicity6/thetaSubavg/ftest_CosSinTheta/spmT_0001.nii')
target_coord = [6, 44, -10];

V = spm_vol(x.fullpath);

voxel_coord = inv(V.mat) * [target_coord 1]';
voxel_coord = round(voxel_coord(1:3));

selected_coord_h = V.mat * [voxel_coord; 1];
selected_coord = selected_coord_h(1:3)';

selected_value = spm_sample_vol(V, voxel_coord(1), voxel_coord(2), voxel_coord(3), 0);

fprintf('Requested MNI coordinate: [%g %g %g]\n', target_coord);
fprintf('Selected MNI coordinate:  [%g %g %g]\n', selected_coord);
fprintf('Voxel indices:             [%g %g %g]\n', voxel_coord);
fprintf('Value: %.4f\n', selected_value);

%Value: 0.1445
%1-tcdf(0.1445,17)
%0.4434
