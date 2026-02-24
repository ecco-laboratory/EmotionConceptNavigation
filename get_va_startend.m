clear; close all; clc;
[project_dir, ~, psychopy_csv_dir, theta_dir, category_dir, vaDiff_dir, avgVA_dir, distance_dir, ~, ~, ~,~, ~, ~, subjects] = set_up_dirs_constants();

va_dir = fullfile(project_dir, 'data', 'beh', 'VA');


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
if ~exist(va_dir, 'dir')
    mkdir(va_dir);
end

for i = 1:length(subjects)
    subject = subjects{i};

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
            % Load norm data
            norm_data_file = fullfile(psychopy_csv_dir, 'behTables',['norm_table_', modality, '.csv']);
            norm_data = readtable(norm_data_file);
            % Load subject-specific ratings
            subject_specific_ratings_file = fullfile(psychopy_csv_dir, 'behTables', ['sub', subject], ['sub', subject, '_table_', modality, '.csv']);
            if exist(subject_specific_ratings_file, 'file')
                subject_specific_ratings = readtable(subject_specific_ratings_file);
            else
                warning(['Could not find subject-specific ratings file for subject ', subject, ', so using sub0001 as template and filling with nans']);
                sub0001_data = readtable(fullfile(psychopy_csv_dir, 'behTables', 'sub0001', ['sub0001_table_', modality, '.csv']));
                subject_specific_ratings = array2table(nan(size(sub0001_data, 1), length(sub0001_data.Properties.VariableNames)), 'VariableNames', sub0001_data.Properties.VariableNames);
            end
            % Load subject-averaged ratings
            subject_averaged_ratings_file = fullfile(psychopy_csv_dir, 'behTables', ['subAvg_table_', modality, '.csv']);
            subject_averaged_ratings = readtable(subject_averaged_ratings_file);
            
            startValence_norm = []; startArousal_norm = [];
            endValence_norm = []; endArousal_norm = [];
            startValence_subspec = []; startArousal_subspec = [];
            endValence_subspec = []; endArousal_subspec = [];
            startValence_subavg = []; startArousal_subavg = [];
            endValence_subavg = []; endArousal_subavg = [];
            
            if strcmp(modality, 'face')
                % For face, find match using face_id, start_category, and end_category
                for trial = 1:size(scanner_beh_file, 1)
                    if isempty(scanner_beh_file.group{trial}) 
                        continue;  % Skip this trial
                    end

                    
                    face_id = scanner_beh_file.group{trial};
                    start_category = scanner_beh_file.start_category{trial};
                    end_category = scanner_beh_file.end_category{trial};
                    %Match in norm data
                    norm_row = find(strcmp(norm_data.face_id, face_id) & ...
                                    strcmp(norm_data.start_label, start_category) & ...
                                    strcmp(norm_data.end_label, end_category));
                    if ~isempty(norm_row)
                        startValence_norm(end+1) = norm_data.start_valence(norm_row);
                        startArousal_norm(end+1) = norm_data.start_arousal(norm_row);
                        endValence_norm(end+1) = norm_data.end_valence(norm_row);
                        endArousal_norm(end+1) = norm_data.end_arousal(norm_row);
                    else
                        startValence_norm(end+1) = NaN;
                        startArousal_norm(end+1) = NaN;
                        endValence_norm(end+1) = NaN;
                        endArousal_norm(end+1) = NaN;
                        warning(['No matching norm data for face_id: ', face_id, ', start_category: ', start_category, ', end_category: ', end_category]);
                    end
                    
                    %Match in subject-specific ratings
                    subspec_row = find(strcmp(subject_specific_ratings.group, face_id) & ...
                                      strcmp(subject_specific_ratings.ori_start_category, start_category) & ...
                                      strcmp(subject_specific_ratings.ori_end_category, end_category));
                    if ~isempty(subspec_row)
                        startValence_subspec(end+1) = subject_specific_ratings.start_valence(subspec_row);
                        startArousal_subspec(end+1) = subject_specific_ratings.start_arousal(subspec_row);
                        endValence_subspec(end+1) = subject_specific_ratings.end_valence(subspec_row);
                        endArousal_subspec(end+1) = subject_specific_ratings.end_arousal(subspec_row);
                    else
                        startValence_subspec(end+1) = NaN;
                        startArousal_subspec(end+1) = NaN;
                        endValence_subspec(end+1) = NaN;
                        endArousal_subspec(end+1) = NaN;
                        warning(['No matching subject-specific data for face_id: ', face_id, ', start_category: ', start_category, ', end_category: ', end_category]);
                    end

                    %Match in subject-averaged ratings
                    subavg_row = find(strcmp(subject_averaged_ratings.group, face_id) & ...
                                     strcmp(subject_averaged_ratings.ori_start_category, start_category) & ...
                                     strcmp(subject_averaged_ratings.ori_end_category, end_category));
                    if ~isempty(subavg_row)
                        startValence_subavg(end+1) = subject_averaged_ratings.start_valence(subavg_row);
                        startArousal_subavg(end+1) = subject_averaged_ratings.start_arousal(subavg_row);
                        endValence_subavg(end+1) = subject_averaged_ratings.end_valence(subavg_row);
                        endArousal_subavg(end+1) = subject_averaged_ratings.end_arousal(subavg_row);
                    else
                        startValence_subavg(end+1) = NaN;
                        startArousal_subavg(end+1) = NaN;
                        endValence_subavg(end+1) = NaN;
                        endArousal_subavg(end+1) = NaN;
                        warning(['No matching subject-averaged data for face_id: ', face_id, ', start_category: ', start_category, ', end_category: ', end_category]);
                    end
                end

                
                if ~isempty(startValence_norm)
                   %check length of startValence_norm, startValence_subspec, startValence_subavg
                   if length(scanner_beh_file.group(2:end-1)) ~= length(startValence_norm) || length(scanner_beh_file.group(2:end-1)) ~= length(startValence_subspec) || length(scanner_beh_file.group(2:end-1)) ~= length(startValence_subavg)
                       warning(['Length of scanner_beh_file(2:end-1) is not the same as the length of thetas_norm, thetas_subspec, thetas_subavg for subject ', subject]);
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
                    
                    % Match in norm data
                    norm_row = find(strcmp(norm_data.group, group) & ...
                                    strcmp(norm_data.start_stimulus, start_word) & ...
                                    strcmp(norm_data.end_stimulus, end_word));
                    if ~isempty(norm_row)
                        startValence_norm(end+1) = norm_data.start_valence(norm_row);
                        endValence_norm(end+1) = norm_data.end_valence(norm_row);
                        startArousal_norm(end+1) = norm_data.start_arousal(norm_row);
                        endArousal_norm(end+1) = norm_data.end_arousal(norm_row);
                    else
                        startValence_norm(end+1) = NaN;
                        endValence_norm(end+1) = NaN;
                        startArousal_norm(end+1) = NaN;
                        endArousal_norm(end+1) = NaN;
                        warning(['No matching norm data for group: ', group, ', start_word: ', start_word, ', end_word: ', end_word]);
                    end
                    
                    % Match in subject-specific ratings
                    subspec_row = find(strcmp(subject_specific_ratings.group, group) & ...
                                      strcmp(subject_specific_ratings.start_stimulus, start_word) & ...
                                      strcmp(subject_specific_ratings.end_stimulus, end_word));
                    if ~isempty(subspec_row)
                        startValence_subspec(end+1) = subject_specific_ratings.start_valence(subspec_row);
                        endValence_subspec(end+1) = subject_specific_ratings.end_valence(subspec_row);
                        startArousal_subspec(end+1) = subject_specific_ratings.start_arousal(subspec_row);
                        endArousal_subspec(end+1) = subject_specific_ratings.end_arousal(subspec_row);

                    else
                        startValence_subspec(end+1) = NaN;
                        endValence_subspec(end+1) = NaN;
                        startArousal_subspec(end+1) = NaN;
                        endArousal_subspec(end+1) = NaN;
                        warning(['No matching subject-specific data for group: ', group, ', start_word: ', start_word, ', end_word: ', end_word]);
                    end
                    
                    % Match in subject-averaged ratings
                    subavg_row = find(strcmp(subject_averaged_ratings.group, group) & ...
                                     strcmp(subject_averaged_ratings.start_stimulus, start_word) & ...
                                     strcmp(subject_averaged_ratings.end_stimulus, end_word));
                    if ~isempty(subavg_row)
                        startValence_subavg(end+1) = subject_averaged_ratings.start_valence(subavg_row);
                        endValence_subavg(end+1) = subject_averaged_ratings.end_valence(subavg_row);
                        startArousal_subavg(end+1) = subject_averaged_ratings.start_arousal(subavg_row);
                        endArousal_subavg(end+1) = subject_averaged_ratings.end_arousal(subavg_row);       
                    else
                        startValence_subavg(end+1) = NaN;
                        endValence_subavg(end+1) = NaN;
                        startArousal_subavg(end+1) = NaN;
                        endArousal_subavg(end+1) = NaN;
                        warning(['No matching subject-averaged data for group: ', group, ', start_word: ', start_word, ', end_word: ', end_word]);
                    end
                end
                

                if ~isempty(startValence_norm)
                   %check length of startValence_norm, startValence_subspec, startValence_subavg
                   if length(scanner_beh_file.group(2:end-1)) ~= length(startValence_norm) || length(scanner_beh_file.group(2:end-1)) ~= length(startValence_subspec) || length(scanner_beh_file.group(2:end-1)) ~= length(startValence_subavg)
                       warning(['Length of scanner_beh_file(2:end-1) is not the same as the length of thetas_norm, thetas_subspec, thetas_subavg for subject ', subject]);
                       continue;
                   end
                end
            end
            

            va = struct();
            va.startValence.Norm = startValence_norm;
            va.startArousal.Norm = startArousal_norm;
            va.endValence.Norm = endValence_norm;
            va.endArousal.Norm = endArousal_norm;
            va.startValence.Subspec = startValence_subspec;
            va.startArousal.Subspec = startArousal_subspec;
            va.endValence.Subspec = endValence_subspec;
            va.endArousal.Subspec = endArousal_subspec;
            va.startValence.Subavg = startValence_subavg;
            va.startArousal.Subavg = startArousal_subavg;
            va.endValence.Subavg = endValence_subavg;
            va.endArousal.Subavg = endArousal_subavg;
            va_file = fullfile(va_dir, ['sub', subject], ['sub', subject, '_', modality, '_va_run', num2str(run), '.mat']);
            if ~exist(fullfile(va_dir, ['sub', subject]), 'dir')
                mkdir(fullfile(va_dir, ['sub', subject]));
            end
            save(va_file, 'va');
            disp(['Saved va to ', va_file]);
        end
    end
end

disp('Processing completed for all subjects, modalities, and runs.');


