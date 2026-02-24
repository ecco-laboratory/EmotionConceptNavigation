[project_dir, fmri_data_dir, ~, theta_dir, ~, ~, ~, ~, spm_timing_dir, motion_regressors_dir, nuisance_regressors_dir, spm1stlevel_dir, ~, ~, subjects] = set_up_dirs_constants();

confounds_to_use = 'motion'; %'motion'; %'motion_rejectedICA'; %'motion'
add_button = 'Button';%'';
modalities = {'word', 'face'};
runs = [1, 2]; % Two runs for each modality
TR = 1.25; % TR in seconds
discard_time = 7.5; % Time to discard in seconds
discard_volumes = discard_time / TR; % Number of volumes to discard 
theta_sources = {'Subavg', 'Subspec'};
periodicity = {6, 4, 5, 7, 8};


for p = 1:length(periodicity)
    current_periodicity = periodicity{p};
    for t = 1:length(theta_sources)
        current_theta_source = theta_sources{t};
        output_dir = fullfile(project_dir, 'outputs', 'bin_trial_counts', ['periodicity', num2str(current_periodicity)], current_theta_source);
        if ~exist(output_dir, 'dir')
            mkdir(output_dir);
        end
        output_file = fullfile(output_dir, 'bin_trial_counts.csv');

        for s = 1:length(subjects)
            subject = subjects{s};
            
            %% LOOP THROUGH MODALITIES AND RUNS TO ADD SESSIONS
            for m = 1:length(modalities)
                modality = modalities{m};
                
                for r = 1:length(runs)
                    run = runs(r);

                    angle_file = fullfile(theta_dir, ['sub', subject], ['sub', subject, '_', modality, '_thetas_run', num2str(run), '.mat']);
                    load(angle_file);
                    angle_data = thetas.(current_theta_source)';
                    angle_avg_data = thetas.Subavg';
                    if strcmp(current_theta_source, 'Subspec')
                        nan_idx = isnan(angle_data);
                        angle_data(nan_idx) = angle_avg_data(nan_idx);
                    end
                    [bin_centers_trials, bin_types_trials] = getGridBinTypes(angle_data, current_periodicity);
                    bin_centers_unique = unique(bin_centers_trials); 

                    for bin_idx = 1:length(bin_centers_unique)
                        current_bin_center = bin_centers_unique(bin_idx); 
                        current_bin_type = unique(bin_types_trials(bin_centers_trials == current_bin_center));
                        if numel(current_bin_type) > 1
                            error('Multiple bin types found for the same bin center');
                        end
                        if current_bin_type == 1, current_bin_type_str = 'on'; elseif current_bin_type == -1, current_bin_type_str = 'off'; else error('Invalid bin type'); end
 
                        current_bin_trial_count = sum(bin_centers_trials == current_bin_center);
                        bin_trial_counts_row = table(current_bin_trial_count, current_bin_center, string(current_bin_type_str), string(modality), run, string(subject), ...
                                                    'VariableNames', {'trial_count', 'bin_center', 'bin_type', 'modality', 'run', 'subject'});
                        if ~exist(output_file, 'file')
                            writetable(bin_trial_counts_row, output_file);
                        else
                            writetable(bin_trial_counts_row, output_file, 'WriteMode', 'append');
                        end
                    end
                    
                    
                end
            end
        end
    end
end


function [bin_centers, bin_types] = getGridBinTypes(angle_data, periodicity)
    %deal with float-point error
    tol = 10 * eps(max(abs(angle_data)));

    % classify angles into: 1 = on-grid, -1 = off-grid
    bin_centers = nan(size(angle_data));
    bin_types = nan(size(angle_data));

    [on_grid_bins, off_grid_bins, on_centers, off_centers] = get_grid_bins_with_edges(periodicity);
    %the first two on_grid_bins are share the same bin_center, so add one more bin_center
    on_centers = [on_centers(1), on_centers];
    
    % assign on-grid = 1
    for i = 1:size(on_grid_bins, 1)
        lower = on_grid_bins(i, 1);
        upper = on_grid_bins(i, 2);
        idx = angle_data >= lower - tol & angle_data < upper + tol & isnan(bin_types);
        bin_types(idx) = 1;
        bin_centers(idx) = on_centers(i);
    end
    
    % assign off-grid = -1
    for i = 1:size(off_grid_bins, 1)
        lower = off_grid_bins(i, 1);
        upper = off_grid_bins(i, 2);
        idx = angle_data >= lower - tol & angle_data < upper + tol & isnan(bin_types);
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
    