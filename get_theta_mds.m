clear; close all; clc;
[project_dir, ~, psychopy_csv_dir, ~, ~, ~, ~, ~, ~, ~, ~,~, ~, ~, subjects] = set_up_dirs_constants();


mds_metrics = {'cosine', 'correlation'};
mds_aligned = true;
if mds_aligned, seeds = {'seedAvg'};else seeds = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10'}; end
if mds_aligned, mds_folder = 'MDS_aligned'; else mds_folder = 'MDS'; end
overwrite = true;
runs = [1, 2];
modalities = {'face', 'word'};
category_names = {'hap','fear','sad','surprise','disgust','neutral','anger'};
cate_ori_rating_names = struct( ...
    'HAP', 'hap', ...
    'AFR', 'fear', ...
    'SAD', 'sad', ...
    'SUR', 'surprise', ...
    'DIS', 'disgust', ...
    'NEU', 'neutral', ...
    'ANG', 'anger');

num_trials_face = 84;
num_trials_word = 124;


for metr = 1:length(mds_metrics)
    metric = mds_metrics{metr};
    if mds_aligned, theta_dir = fullfile(project_dir, 'data', 'beh', ['theta_mdsAligned_', metric]); else theta_dir = fullfile(project_dir, 'data', 'beh', ['theta_mds_', metric]); end
    if ~exist(theta_dir, 'dir') mkdir(theta_dir); end
    for see = 1:length(seeds)
        seed = seeds{see};
        if isnumeric(seed), seed = num2str(seed); end
        for i = 1:length(subjects)
            subject = subjects{i};
            %check if the file already exists
            if ~overwrite && exist(fullfile(theta_dir, ['seed', seed], ['sub', subject], ['sub', subject, '_face_thetas_run1.mat']), 'file') && ...
            exist(fullfile(theta_dir, ['seed', seed], ['sub', subject], ['sub', subject, '_face_thetas_run2.mat']), 'file') && ...
            exist(fullfile(theta_dir, ['seed', seed], ['sub', subject], ['sub', subject, '_word_thetas_run1.mat']), 'file') && ...
            exist(fullfile(theta_dir, ['seed', seed], ['sub', subject], ['sub', subject, '_word_thetas_run2.mat']), 'file') 
            
                disp(['Skipping subject ', subject, ' because files already exist']);
                continue;
            end
            for m = 1:length(modalities)
                modality = modalities{m};
                for r = 1:length(runs)
                    run = runs(r);
                    disp(['Processing subject: ', subject, ', run: ', num2str(run), ', modality: ', modality]);
                    
                    scanner_beh_files = dir(fullfile(psychopy_csv_dir, ['sub-', subject], ['sub-', subject, '_task-', modality, '_run-0', num2str(run), '_*.csv']));
                    if isempty(scanner_beh_files) 
                        warning(['Could not find scanner behavioral files for subject ', subject]);
                        continue;
                    end
                    scanner_beh_file = readtable(fullfile(scanner_beh_files(1).folder, scanner_beh_files(1).name));
                    % Load subject-specific ratings
                    subject_specific_ratings_file = fullfile(psychopy_csv_dir, 'behTables', mds_folder, ['sub', subject], ['sub', subject, '_table_', modality, '_mds', metric, seed, '.csv']);
                    if exist(subject_specific_ratings_file, 'file')
                        subject_specific_ratings = readtable(subject_specific_ratings_file);
                    else
                        warning(['Could not find subject-specific ratings file for subject ', subject, ', so using sub0001 as template and filling with nans']);
                        sub0001_data = readtable(fullfile(psychopy_csv_dir, 'behTables', mds_folder, 'sub0001', ['sub0001_table_', modality, '_mds', metric, seed, '.csv']));
                        subject_specific_ratings = array2table(nan(size(sub0001_data, 1), length(sub0001_data.Properties.VariableNames)), 'VariableNames', sub0001_data.Properties.VariableNames);
                    end
                    % Load subject-averaged ratings
                    subject_averaged_ratings_file = fullfile(psychopy_csv_dir, 'behTables', mds_folder, ['subAvg_table_', modality, '_mds', metric, seed, '.csv']);
                    subject_averaged_ratings = readtable(subject_averaged_ratings_file);
                    
                
                    thetas_subspec = [];
                    thetas_subavg = [];
                    
                    if strcmp(modality, 'face')
                        % For face, find match using face_id, start_category, and end_category
                        for trial = 1:size(scanner_beh_file, 1)
                            if isempty(scanner_beh_file.group{trial}) 
                                continue;  % Skip this trial
                            end

                            
                            face_id = scanner_beh_file.group{trial};
                            start_category = scanner_beh_file.start_category{trial};
                            end_category = scanner_beh_file.end_category{trial};
                            
                            %Match in subject-specific ratings
                            subspec_row = find(strcmp(subject_specific_ratings.group, face_id) & ...
                                            strcmp(subject_specific_ratings.ori_start_category, start_category) & ...
                                            strcmp(subject_specific_ratings.ori_end_category, end_category));
                            if ~isempty(subspec_row)
                                thetas_subspec(end+1) = subject_specific_ratings.theta(subspec_row);
                            else
                                thetas_subspec(end+1) = NaN;
                            end

                            %Match in subject-averaged ratings
                            subavg_row = find(strcmp(subject_averaged_ratings.group, face_id) & ...
                                            strcmp(subject_averaged_ratings.ori_start_category, start_category) & ...
                                            strcmp(subject_averaged_ratings.ori_end_category, end_category));
                            if ~isempty(subavg_row)
                                thetas_subavg(end+1) = subject_averaged_ratings.theta(subavg_row);
                            else
                                thetas_subavg(end+1) = NaN;
                            end
                        end

                        
                        if ~isempty(thetas_subspec)
                        %check length of thetas_subspec, thetas_subavg
                        if length(scanner_beh_file.group(2:end-1)) ~= length(thetas_subspec) || length(scanner_beh_file.group(2:end-1)) ~= length(thetas_subavg)
                            warning(['Length of scanner_beh_file(2:end-1) is not the same as the length of thetas_subspec, thetas_subavg for subject ', subject]);
                            continue;
                        end
                        end
                        
                    elseif strcmp(modality, 'word')
                        % For word, find match using group, start_word, and end_word
                        for trial = 1:size(scanner_beh_file, 1)
                            if isempty(scanner_beh_file.group{trial})
                                continue;  % Skip this trial
                            end
                            
                            group = scanner_beh_file.group{trial};
                            start_word = scanner_beh_file.start_word{trial};
                            end_word = scanner_beh_file.end_word{trial};
                            start_category = scanner_beh_file.start_category{trial};
                            end_category = scanner_beh_file.end_category{trial};
                            if isempty(start_word) || isempty(end_word)
                                continue;
                            end
                            
                            % Match in subject-specific ratings
                            subspec_row = find(strcmp(subject_specific_ratings.group, group) & ...
                                            strcmp(subject_specific_ratings.start_stimulus, start_word) & ...
                                            strcmp(subject_specific_ratings.end_stimulus, end_word));
                            if ~isempty(subspec_row)
                                thetas_subspec(end+1) = subject_specific_ratings.theta(subspec_row);
                            
                            else
                                thetas_subspec(end+1) = NaN;
                            end
                            
                            % Match in subject-averaged ratings
                            subavg_row = find(strcmp(subject_averaged_ratings.group, group) & ...
                                            strcmp(subject_averaged_ratings.start_stimulus, start_word) & ...
                                            strcmp(subject_averaged_ratings.end_stimulus, end_word));
                            if ~isempty(subavg_row)
                                thetas_subavg(end+1) = subject_averaged_ratings.theta(subavg_row);
                                
                            else
                                thetas_subavg(end+1) = NaN;
                            
                            end
                        end
                        

                        if ~isempty(thetas_subspec)
                        %check length of thetas_subspec, thetas_subavg
                        if length(scanner_beh_file.group(2:end-1)) ~= length(thetas_subspec) || length(scanner_beh_file.group(2:end-1)) ~= length(thetas_subavg)
                            warning(['Length of scanner_beh_file(2:end-1) is not the same as the length of thetas_subspec, thetas_subavg for subject ', subject]);
                            continue;
                        end
                        end
                    end
                    
                    thetas = struct();
                    thetas.Subspec = thetas_subspec;
                    thetas.Subavg = thetas_subavg;
                    theta_file = fullfile(theta_dir, ['seed', seed], ['sub', subject], ['sub', subject, '_', modality, '_thetas_run', num2str(run), '.mat']);
                    if ~exist(fullfile(theta_dir, ['seed', seed], ['sub', subject]), 'dir')
                        mkdir(fullfile(theta_dir, ['seed', seed], ['sub', subject]));
                    end
                    save(theta_file, 'thetas');
                    disp(['Saved thetas to ', theta_file]);
                end
            end
        end
    end
end
disp('Processing completed for all subjects, modalities, and runs.');
