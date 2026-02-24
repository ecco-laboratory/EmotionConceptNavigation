[project_dir, fmri_data_dir, ~, theta_dir, ~, ~, ~, ~, spm_timing_dir, motion_regressors_dir, nuisance_regressors_dir, spm1stlevel_dir, ~, phi_dir, subjects] = set_up_dirs_constants();

ICA_type = 'noICA'; 
spm1stlevel_dir = [spm1stlevel_dir, '_', ICA_type, '_wButton'];
theta_sources = {'Subavg'};%, 'Subspec'};
periodicity = {6};%, 4, 5, 7, 8};
cross_validations = {'xModalityRun', 'xModality'};%, 'xRun'};
phi_averaging_methods = 'voxelComponentAverage'; 
phi_calculation_methods = 'singleTrialBeta'; 
phi_source_region = 'current';
base_output_dir = fullfile(project_dir, 'outputs', 'singleTrialBetaAnalysis', ICA_type, 'incl_all_subs_trials','onoffGridcontrast_multivariate_searchlight',['phi_',phi_source_region]);
if ~exist(base_output_dir, 'dir'), mkdir(base_output_dir); end
runs = [1, 2]; 
searchlight_radius = 5;


brain_atlas = load_atlas('canlab2018');
region_names = {'wholebrain'};
region_masks = {[]};

template_nifti_object = fmri_data(fullfile(spm1stlevel_dir, ['sub', subjects{1}], 'singleTrial', 'beta_0001.nii'));
[brain_data_wholebrain, brain_modality_run_idx, brain_sub_ids] = load_and_prepare_brain_data(subjects, spm1stlevel_dir);


for t = 1:length(theta_sources)
    current_angle_source = theta_sources{t};
    [angle_data_all, subject_ids, modality_run_ids] = load_and_prepare_angle_data(subjects, runs, spm1stlevel_dir, theta_dir, current_angle_source);

    for p = 1:length(periodicity)
        current_periodicity = periodicity{p};
        fprintf('Processing periodicity %d...\n', current_periodicity);
        for s = 1:length(subjects)
            current_subject = subjects{s}; current_subject_idx = strcmp(subject_ids, current_subject);
    
            for cv = 1:length(cross_validations)
                current_cross_validation = cross_validations{cv};

                if strcmp(current_cross_validation, 'xModalityRun')
                    training_idx = {{'face1', 'word1'}, {'face1', 'word2'}, {'face2', 'word1'}, {'face2', 'word2'}};
                    test_idx = {{'face2', 'word2'}, {'face2', 'word1'}, {'face1', 'word2'}, {'face1', 'word1'}};
                elseif strcmp(current_cross_validation, 'xModality')
                    training_idx = {{'face1','face2'}, {'word1', 'word2'}};
                    test_idx = {{'word1', 'word2'}, {'face1', 'face2'}};
                elseif strcmp(current_cross_validation, 'xRun')
                    training_idx = {{'face1'}, {'face2'}, {'word1'}, {'word2'}};
                    test_idx = {{'face2'}, {'face1'}, {'word2'}, {'word1'}};
                end

                for t_idx = 1:length(training_idx)
                    current_training_idx = training_idx{t_idx}; current_test_idx = test_idx{t_idx};
                
                    for reg = 1:length(region_names)
                        if strcmp(region_names{reg}, 'wholebrain')
                            searchlight_data = template_nifti_object;
                        else
                            searchlight_data = apply_mask(template_nifti_object, region_masks{reg});
                        end
                        results_vec = run_searchlight(searchlight_data, @perform_multivariate_contrast, 'r', searchlight_radius, ...
                                            'singleTrial_betas', brain_data_wholebrain(current_subject_idx & ismember(modality_run_ids, current_test_idx), :), ...
                                            'angle_data', angle_data_all(current_subject_idx & ismember(modality_run_ids, current_test_idx), :), ...
                                            'singleTrial_betas_phi', brain_data_wholebrain(current_subject_idx & ismember(modality_run_ids, current_training_idx), :), ...
                                            'angle_data_phi', angle_data_all(current_subject_idx & ismember(modality_run_ids, current_training_idx), :), ...
                                            'current_periodicity', current_periodicity);

                        %[contrast_value, bin_means, bin_centers, num_aligned_trials] = perform_multivariate_contrast(brain_data_wholebrain(current_subject_idx & ismember(modality_run_ids, current_test_idx), :), ...
                        %                                           angle_data_all(current_subject_idx & ismember(modality_run_ids, current_test_idx), :), phi_value, current_periodicity);
                        path_to_save = fullfile(base_output_dir, current_angle_source, ['periodicity', num2str(current_periodicity)], current_cross_validation, ['phi_', phi_averaging_methods], phi_calculation_methods, 'nifti');
                        if ~exist(path_to_save, 'dir'), mkdir(path_to_save); end
                        searchlight_data.dat = results_vec;
                        searchlight_data.fullpath = fullfile(path_to_save, ['contrast_value_sub', current_subject, '_', strjoin(current_test_idx, '_'), '_', region_names{reg}, '.nii']);
                        searchlight_data.write;
                    end
                end
            end   
        end
    end
end

function [brain_data_wholebrain, brain_modality_run_idx, brain_sub_ids] = load_and_prepare_brain_data(subjects, spm1stlevel_dir)
    beta_images = {};
    brain_modality_run_idx = {};
    brain_sub_ids = {};

    for s = 1:length(subjects)
        current_subject = subjects{s};
        
        singleTrial_dir = fullfile(spm1stlevel_dir, ['sub', current_subject], 'singleTrial');
        load(fullfile(singleTrial_dir, 'SPM.mat'));
        spm_regressor_names = SPM.xX.name;
        stim_trials = contains(spm_regressor_names, 'trial');
        stim_trials_idx = find(stim_trials);
        
        % Prepare beta images
        beta_images_current_subject = {};
        brain_modality_run_idx_current_subject = {};
        brain_sub_ids_current_subject = {};
        for i = 1:length(stim_trials_idx)
            current_trial = stim_trials_idx(i);
            current_beta_image = fullfile(singleTrial_dir, ['beta_', num2str(current_trial, '%04d'), '.nii']);
            beta_images_current_subject{i, 1} = current_beta_image;
            stim_name = spm_regressor_names{current_trial};
            stim_info = regexp(stim_name, 'Sn\((\d+)\)\s+(\w+)_run(\d+)_trial(\d+)\*bf\((\d+)\)', 'tokens');
            brain_modality_run_idx_current_subject{i, 1} = [stim_info{1}{2}, stim_info{1}{3}];
            brain_sub_ids_current_subject{i, 1} = current_subject;
        end
        beta_images = [beta_images; beta_images_current_subject];
        brain_modality_run_idx = [brain_modality_run_idx; brain_modality_run_idx_current_subject];
        brain_sub_ids = [brain_sub_ids; brain_sub_ids_current_subject];
    end
    % Load brain data
    brain_data_wholebrain = fmri_data(beta_images).dat';
end

function [angle_data_all, subject_ids, modality_run_ids] = load_and_prepare_angle_data(subjects, runs, spm1stlevel_dir, theta_dir, current_angle_source)

    angle_data_all = [];
    subject_ids = {};
    modality_run_ids = {};

    for s = 1:length(subjects)
        current_subject = subjects{s};
        
        singleTrial_dir = fullfile(spm1stlevel_dir, ['sub', current_subject], 'singleTrial');
        load(fullfile(singleTrial_dir, 'SPM.mat'));
        spm_regressor_names = SPM.xX.name;
        stim_trials = contains(spm_regressor_names, 'trial');
        stim_trials_idx = find(stim_trials);
        
        first_stim_trial = spm_regressor_names(find(stim_trials, 1, 'first'));
        last_stim_trial = spm_regressor_names(find(stim_trials, 1, 'last'));
        
        if contains(first_stim_trial, 'face') && contains(last_stim_trial, 'word')
            modalities = {'face', 'word'};
        elseif contains(first_stim_trial, 'word') && contains(last_stim_trial, 'face')
            modalities = {'word', 'face'};
        else
            error('Modality order unclear.');
        end
        
        % Load angle data
        angles_all_modalities_runs = [];
        angles_modality_run_idx = {};
        for m = 1:length(modalities)
            current_modality = modalities{m};
            for r = 1:length(runs)
                current_run = runs(r);
                angle_file = fullfile(theta_dir, ['sub', current_subject], ['sub', current_subject, '_', current_modality, '_thetas_run', num2str(current_run), '.mat']);
                load(angle_file);
                
                angle_data = thetas.(current_angle_source)';
                angle_avg_data = thetas.Subavg';
                if strcmp(current_angle_source, 'Subspec')
                    nan_idx = isnan(angle_data);
                    angle_data(nan_idx) = angle_avg_data(nan_idx);
                end
                
                angles_all_modalities_runs = [angles_all_modalities_runs; angle_data];
                current_modality_run_idx = [current_modality, num2str(current_run)];
                angles_modality_run_idx = [angles_modality_run_idx; repmat({current_modality_run_idx}, size(angle_data, 1), 1)];
            end
        end
        angle_data_all = [angle_data_all; angles_all_modalities_runs];
        subject_ids = [subject_ids; repmat({current_subject}, size(angles_all_modalities_runs, 1), 1)];
        modality_run_ids = [modality_run_ids; angles_modality_run_idx];
    end
end

function [bin_centers, bin_types] = getGridBinTypes(aligned_angle_data, periodicity)
    %deal with float-point error
    tol = 10 * eps(max(abs(aligned_angle_data)));

    % classify angles into: 1 = on-grid, -1 = off-grid
    bin_centers = nan(size(aligned_angle_data));
    bin_types = nan(size(aligned_angle_data));

    [on_grid_bins, off_grid_bins, on_centers, off_centers] = get_grid_bins_with_edges(periodicity);
    %the first two on_grid_bins are share the same bin_center, so add one more bin_center
    on_centers = [on_centers(1), on_centers];
    
    % assign on-grid = 1
    for i = 1:size(on_grid_bins, 1)
        lower = on_grid_bins(i, 1);
        upper = on_grid_bins(i, 2);
        idx = aligned_angle_data >= lower - tol & aligned_angle_data < upper + tol & isnan(bin_types);
        bin_types(idx) = 1;
        bin_centers(idx) = on_centers(i);
    end
    
    % assign off-grid = -1
    for i = 1:size(off_grid_bins, 1)
        lower = off_grid_bins(i, 1);
        upper = off_grid_bins(i, 2);
        idx = aligned_angle_data >= lower - tol & aligned_angle_data < upper + tol & isnan(bin_types);
        bin_types(idx) = -1;
        bin_centers(idx) = off_centers(i);
    end
    
end
   

function [on_grid_bins, off_grid_bins, on_centers, off_centers] = get_grid_bins_with_edges(periodicity)
    
        bin_width = 360 / (2 * periodicity);                    
        num_bins = 2 * periodicity;
    
        grid_angles = mod((0:periodicity-1) * 360/periodicity, 360);  
        on_centers = grid_angles;                                        
        off_centers = mod(grid_angles + 360/(2 * periodicity), 360);           
    
        on_grid_bins = expand_bins(on_centers, bin_width);
        off_grid_bins = expand_bins(off_centers, bin_width);
end
    
function bins = expand_bins(centers, width)
    
    half_width = width / 2;
    bins = [];

    for i = 1:length(centers)
        start_angle = mod(centers(i) - half_width, 360);
        end_angle = mod(centers(i) + half_width, 360);

        if start_angle < end_angle
            bins = [bins; start_angle, end_angle];
        else
            bins = [bins; start_angle, 360; 0, end_angle];
        end
    end
end

function phi_value = cal_phi_value(singleTrial_betas_phi, angle_data_phi, current_periodicity)
    beta_sin_allvoxels = nan(size(singleTrial_betas_phi, 2), 1);
    beta_cos_allvoxels = nan(size(singleTrial_betas_phi, 2), 1);
    X_train = [cos(current_periodicity * deg2rad(angle_data_phi)), sin(current_periodicity * deg2rad(angle_data_phi))];
    X_design = [ones(size(X_train, 1), 1), X_train];
    for v = 1:size(singleTrial_betas_phi, 2)
        y = singleTrial_betas_phi(:, v);
        b_noPhi = X_design \ y;
        beta_cos_allvoxels(v) = b_noPhi(2);
        beta_sin_allvoxels(v) = b_noPhi(3);
    end
    phi_value = atan2(mean(beta_sin_allvoxels), mean(beta_cos_allvoxels)) / current_periodicity;
end

function contrast_value = perform_multivariate_contrast(singleTrial_betas, angle_data, singleTrial_betas_phi,angle_data_phi,current_periodicity)
        phi_value = cal_phi_value(singleTrial_betas_phi, angle_data_phi, current_periodicity);

        phi_deg = rad2deg(phi_value); % convert phi to degrees

        %define bin width and centers
        period_deg = 360 / current_periodicity; %60 for 6-fold
        bin_width = period_deg/2; %30 for 6-fold
        bin_half = bin_width/2; %15 for 6-fold
        bin_centers = 0:bin_width:(180-bin_width); 
        align_half = bin_half/2; %7.5 for 6-fold

        %find aligned trials
        angle_modP = mod(angle_data - phi_deg, period_deg);
        aligned_trials = (abs(angle_modP) < align_half) | (abs(angle_modP - period_deg) < align_half); %deal with wrap-around (like 59-degree should belong to 0 degree bin for 6 fold)
        aligned_idx = find(aligned_trials); 

        nTrials = size(singleTrial_betas, 1);

        % compute correlations for aligned × others
        aligned_corrs = [];
        misaligned_corrs = [];
        bin_corrs = cell(1, length(bin_centers));
        for b = 1:length(bin_centers), bin_corrs{b} = []; end

        for i = 1:numel(aligned_idx)
            t = aligned_idx(i); 
            beta_t = singleTrial_betas(t,:)'; %column vector for trial t

            for j = 1:nTrials
                if j == t, continue; end %skip self

                r = corr(beta_t, singleTrial_betas(j,:)');

                % angular difference modulo 60
                d_raw = abs(angle_data(t) - angle_data(j)); d_on_period = mod(d_raw, period_deg);
                diff_to_0 = min(d_on_period, period_deg - d_on_period);
                diff_to_half = abs(d_on_period - period_deg/2); diff_to_half = min(diff_to_half, period_deg - diff_to_half);
                % classify correlation
                if diff_to_0 < bin_half
                    aligned_corrs(end+1) = r;
                elseif diff_to_half < bin_half
                    misaligned_corrs(end+1) = r;
                end

                %6-bin similarity
                d180 = mod(d_raw, 180);
                for b = 1:length(bin_centers)
                    diff_bincenter = abs(d180 - bin_centers(b)); diff_bincenter = min(diff_bincenter, 180 - diff_bincenter);
                    if diff_bincenter < bin_half
                        bin_corrs{b} = [bin_corrs{b}; r];
                        break;
                    end
                end


            end
        end

        %contrast (aligned vs misaligned)
        contrast_value = mean(aligned_corrs,'omitnan') - mean(misaligned_corrs,'omitnan');

end



function results_vec = run_searchlight(searchlight_data, custom_function, varargin)
    % Run searchlight analysis with a custom function
    %
    % Usage:
    %   results_vec = run_searchlight(searchlight_data, @perform_multivariate_contrast, ...
    %                'r', 5, 'angle_data', angle_data, 'phi_value', phi, 'current_periodicity', 6);
    %
    % Inputs:
    %   searchlight_data : fmri_data object
    %
    %   custom_function : function handle
    %       A function of the form:
    %           outval = custom_function(singleTrial_betas, varargin{:})
    %
    % Optional Inputs:
    %   'r'        : searchlight radius in mm (default = 3)
    %   'indx'     : (optional) precomputed sphere index matrix (nVox × nSpheres, logical)
    %   Any other inputs will be passed through to custom_function
    %
    % Output:
    %   results_vec : vector (nVoxels × 1)
    %       Each voxel’s value is the custom_function result for its sphere.
    %
    
    % -------------------------
    % Parse inputs
    % -------------------------
    r = [];
    indx = [];
    singleTrial_betas = []; angle_data = []; singleTrial_betas_phi = []; angle_data_phi = []; current_periodicity = [];
    xyz = searchlight_data.volInfo.xyzlist;
    %nvox = searchlight_data.volInfo.n_inmask;
    v_size = abs(diag(searchlight_data.volInfo.mat(1:3,1:3)))';
    if isempty(searchlight_data.removed_voxels) || isscalar(searchlight_data.removed_voxels)
        voxs_to_keep = 1:size(searchlight_data.dat,1);
     else
        voxs_to_keep = find(~searchlight_data.removed_voxels);
     end

    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                case 'r'
                    r = varargin{i+1};
                case 'indx'
                    indx = varargin{i+1};
                case 'singleTrial_betas'
                    singleTrial_betas = varargin{i+1};
                case 'angle_data'
                    angle_data = varargin{i+1};
                case 'singleTrial_betas_phi'
                    singleTrial_betas_phi = varargin{i+1};
                case 'angle_data_phi'
                    angle_data_phi = varargin{i+1};
                case 'current_periodicity'
                    current_periodicity = varargin{i+1};
            end
        end
    end
    
    
    % -------------------------
    % Build sphere indices
    % -------------------------
    %if isempty(indx)
        %indx = searchlight_sphere_prep(searchlight_data, r);
    %else
        %fprintf('Using provided searchlight spheres...\n');
    %end
    
    % -------------------------
    % Run full searchlight
    % -------------------------
    fprintf('Running searchlight with %d voxels...\n', length(voxs_to_keep));
    results_vec = nan(length(voxs_to_keep), 1);
    
    tic;
    update_every = 10000;
    for ii = 1:length(voxs_to_keep)
        v = voxs_to_keep(ii);
        center = xyz(v, :);
        %sphere_vox = sum((xyz - center).^2, 2) <= r^2; %this is r in voxels
        sphere_vox = sum(((xyz - center) .* v_size).^2, 2) <= r^2; %this is r in mm
        if ~any(sphere_vox), continue; end
        sphere_data = singleTrial_betas(:, sphere_vox);
        sphere_data_phi = singleTrial_betas_phi(:, sphere_vox);
        try
            results_vec(ii) = custom_function(sphere_data, angle_data, sphere_data_phi, angle_data_phi, current_periodicity);
        catch
            results_vec(ii) = NaN;
            fprintf('Error in custom_function for voxel %d\n', v);
        end
        if mod(ii, update_every)==0 || ii==length(voxs_to_keep)
            elapsed = toc;
            pct = 100 * ii / length(voxs_to_keep);
            avg_time = elapsed / ii;
            remaining = avg_time * (length(voxs_to_keep) - ii);
            fprintf('%.1f%% complete | elapsed %.1f min | remaining %.1f min\n', pct, elapsed/60, remaining/60);
        end
    end
    elapsed = toc;
    
    [hour, minute, second] = sec2hms(elapsed);
    fprintf('Done in %d hr %d min %.1f sec\n', hour, minute, second);
    
end
    
% -------------------------
% Helper: searchlight spheres
% -------------------------
function indx = searchlight_sphere_prep(dat, r)
    
    nvox = dat.volInfo.n_inmask;
    indx = cell(1, nvox);
    
    t=tic;
    fprintf('Preparing %d seeds...\n', nvox);
    
    parfor i = 1:nvox
        seed{i} = dat.volInfo.xyzlist(i, :);
    end
    e = toc(t);
    fprintf('Done in %3.2f sec\n', e);

    % Set up indices for spherical searchlight
    % -------------------------------------------------------------------------
    % These could be indices for ROIs, user input, previously saved indices...

    % First, a rough time estimate:
    % -------------------------------------------------------------------------
    fprintf('Searchlight sphere construction can take 20 mins or more! (est: 20 mins with 8 processors/gray matter mask)\n');
    fprintf('It can be re-used once created for multiple analyses with the same region definitions\n');
    fprintf('Getting a rough time estimate for how long this will take...\n');

    n_to_run = min(500, nvox);
    t = tic;
    parfor i = 1:n_to_run
        
        mydist = sum([dat.volInfo.xyzlist(:, 1) - seed{i}(1) dat.volInfo.xyzlist(:, 2) - seed{i}(2) dat.volInfo.xyzlist(:, 3) - seed{i}(3)] .^ 2, 2);
        indx{i} = mydist <= r.^2;
        
    end
    e = toc(t);
    estim = e * nvox / n_to_run;

    [hour, minute, second] = sec2hms(estim);
    fprintf(1,'\nEstimate for whole brain = %3.0f hours %3.0f min %2.0f sec\n',hour, minute, second);

    % Second, do it for all voxels/spheres:
    % -------------------------------------------------------------------------
    t = tic;
    fprintf('Constructing spheres...\n');
    parfor i = 1:nvox
        mydist = sum([dat.volInfo.xyzlist(:,1)-seed{i}(1), ...
                      dat.volInfo.xyzlist(:,2)-seed{i}(2), ...
                      dat.volInfo.xyzlist(:,3)-seed{i}(3)].^2, 2);
        indx{i} = mydist <= r^2;
    end
    e = toc(t);
    fprintf('Done in %3.2f sec\n', e);

    %Make sparse matrix
    t = tic;
    fprintf('Making sparse matrix...\n');
    indx = sparse(cat(2, indx{:}));

    e = toc(t);
    fprintf('Done in %3.2f sec\n', e);
    
end %end of function searchlight_sphere_prep
    
% -------------------------
% Helper: seconds → h/m/s
% -------------------------
function [hour, minute, second] = sec2hms(sec)
    hour   = fix(sec/3600);
    sec    = sec - 3600*hour;
    minute = fix(sec/60);
    sec    = sec - 60*minute;
    second = sec;
end
 