clear; close all; clc;
[~, ~, psychopy_csv_dir, ~, ~, ~, ~, ~, spm_timing_dir, ~, ~, ~, ~, ~, subjects] = set_up_dirs_constants();

add_button = true;
output_dir = spm_timing_dir;
discard_time = 7.5; % Time to discard in seconds in the beginning of the run
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

for s = 1:length(subjects)
    subject = subjects{s};
    disp(['Processing subject: ', subject]);
    %check if the file already exists
    if exist(fullfile(output_dir, ['sub', subject], ['sub', subject, '_wordRatingButton_run1.mat']), 'file') && ...
       exist(fullfile(output_dir, ['sub', subject], ['sub', subject, '_wordRatingButton_run2.mat']), 'file')
        disp(['Skipping subject ', subject, ' because it is already processed']);
        continue;
    end
    
    %% MAKE SURE THE BEHAVIORAL DATA IS CORRECT
    % Read the csv files (using dir to handle wildcards)
    run1_files = dir(fullfile(psychopy_csv_dir, ['sub-', subject], ['sub-', subject, '*_task-word_run-01_*.csv']));
    run2_files = dir(fullfile(psychopy_csv_dir, ['sub-', subject], ['sub-', subject, '*_task-word_run-02_*.csv']));
    
    if isempty(run1_files) || isempty(run2_files)
        warning(['Could not find behavioral files for subject ', subject]);
        continue;
    end
    psychopy_csv_file_run1 = readtable(fullfile(run1_files(1).folder, run1_files(1).name));
    psychopy_csv_file_run2 = readtable(fullfile(run2_files(1).folder, run2_files(1).name));
    % Count valid trials
    num_trials_run1 = sum(~ismissing(psychopy_csv_file_run1.start_word) & ~isempty(psychopy_csv_file_run1.start_word));
    num_trials_run2 = sum(~ismissing(psychopy_csv_file_run2.start_word) & ~isempty(psychopy_csv_file_run2.start_word));
    % Make sure total trials is 124
    if num_trials_run1 + num_trials_run2 ~= 124
        warning(['Total number of trials is not 124 for words for subject ', subject]);
        % Continue instead of error to process other subjects
        continue;
    end
    % Make sure start_word, iti_cross_started in the first and last row is empty
    if ~ismissing(psychopy_csv_file_run1.start_word(end)) || ~ismissing(psychopy_csv_file_run1.start_word(1)) || ...
       ~ismissing(psychopy_csv_file_run1.iti_cross_started(1)) || ~ismissing(psychopy_csv_file_run1.iti_cross_started(end))
        warning(['start_word, iti_cross_started last and first rows are not empty for run1 for subject ', subject]);
        continue;
    end
    if ~ismissing(psychopy_csv_file_run2.start_word(end)) || ~ismissing(psychopy_csv_file_run2.start_word(1)) || ... 
       ~ismissing(psychopy_csv_file_run2.iti_cross_started(1)) || ~ismissing(psychopy_csv_file_run2.iti_cross_started(end)) 
        warning(['start_word, iti_cross_started last and first rows are not empty for run2 for subject ', subject]);
        continue;
    end
    
    % Get the timing of the scan start
    scan_start_run1 = psychopy_csv_file_run1.trigger_started(1) + psychopy_csv_file_run1.trigger_rt(1) + discard_time; %discard first 8 seconds of the run
    scan_start_run2 = psychopy_csv_file_run2.trigger_started(1) + psychopy_csv_file_run2.trigger_rt(1) + discard_time; %discard first 8 seconds of the run

    %% GET TIMING OF STIMULI
    % Get timing of stimuli
    stim_onsets_run1 = psychopy_csv_file_run1.start_word_text_started(2:end-1) - scan_start_run1;
    stim_onsets_run2 = psychopy_csv_file_run2.start_word_text_started(2:end-1) - scan_start_run2;
    % Move cross_final_started in the last row to the iti_cross_started in the last row
    psychopy_csv_file_run1.iti_cross_started(end) = psychopy_csv_file_run1.cross_final_started(end);
    psychopy_csv_file_run2.iti_cross_started(end) = psychopy_csv_file_run2.cross_final_started(end);
    % Get the duration of the morph video
    stim_durations_run1 = psychopy_csv_file_run1.iti_cross_started(3:end) - psychopy_csv_file_run1.start_word_text_started(2:end-1);
    stim_durations_run2 = psychopy_csv_file_run2.iti_cross_started(3:end) - psychopy_csv_file_run2.start_word_text_started(2:end-1);
    
    %% GET TIMING OF RATING 
    rating_onsets_run1 = psychopy_csv_file_run1.comparison_q_started(3:end-1) - scan_start_run1;
    rating_onsets_run2 = psychopy_csv_file_run2.comparison_q_started(3:end-1) - scan_start_run2;
    rating_durations_run1 = psychopy_csv_file_run1.comparison_q_stopped(3:end-1) - psychopy_csv_file_run1.comparison_q_started(3:end-1);
    rating_durations_run2 = psychopy_csv_file_run2.comparison_q_stopped(3:end-1) - psychopy_csv_file_run2.comparison_q_started(3:end-1);
    % Add the final rating
    rating_onsets_run1 = [rating_onsets_run1; psychopy_csv_file_run1.comparison_q_final_started(end)-scan_start_run1];
    rating_onsets_run2 = [rating_onsets_run2; psychopy_csv_file_run2.comparison_q_final_started(end)-scan_start_run2];
    rating_durations_run1 = [rating_durations_run1; psychopy_csv_file_run1.end_text_started(end)-psychopy_csv_file_run1.comparison_q_final_started(end)];
    rating_durations_run2 = [rating_durations_run2; psychopy_csv_file_run2.end_text_started(end)-psychopy_csv_file_run2.comparison_q_final_started(end)];

    %% GET TIMING OF CONTEXT WORD
    context_onsets_run1 = [];
    context_durations_run1 = [];
    context_onsets_run2 = [];
    context_durations_run2 = [];
    %iterate each row and if it's not empty in the context_text_started column, then add the time to the onsets
    row_idx = 1;
    for i = 1:size(psychopy_csv_file_run1, 1)
        %if the context_text_started is not nan
        if ~ismissing(psychopy_csv_file_run1.context_text_started(i)) 
            context_onsets_run1(row_idx) = psychopy_csv_file_run1.context_text_started(i) - scan_start_run1;
            context_durations_run1(row_idx) = psychopy_csv_file_run1.context_text_stopped(i) - psychopy_csv_file_run1.context_text_started(i);
            row_idx = row_idx + 1;
        end
    end
    row_idx = 1;
    for i = 1:size(psychopy_csv_file_run2, 1)
        if ~ismissing(psychopy_csv_file_run2.context_text_started(i))
            context_onsets_run2(row_idx) = psychopy_csv_file_run2.context_text_started(i) - scan_start_run2;
            context_durations_run2(row_idx) = psychopy_csv_file_run2.context_text_stopped(i) - psychopy_csv_file_run2.context_text_started(i);
            row_idx = row_idx + 1;
        end
    end
    
    % Verify dimensions
    if length(rating_onsets_run1) ~= length(stim_onsets_run1) || length(rating_onsets_run2) ~= length(stim_onsets_run2) || ...
       length(rating_durations_run1) ~= length(stim_durations_run1) || length(rating_durations_run2) ~= length(stim_durations_run2)
        warning(['Length of rating and stim are not the same for subject ', subject]);
        continue;
    end
    names = {'word', 'rating', 'context'};
    onsets_run1 = {stim_onsets_run1, rating_onsets_run1, context_onsets_run1};
    durations_run1 = {stim_durations_run1, rating_durations_run1, context_durations_run1};
    onsets_run2 = {stim_onsets_run2, rating_onsets_run2, context_onsets_run2};
    durations_run2 = {stim_durations_run2, rating_durations_run2, context_durations_run2};
    if ~exist(fullfile(output_dir, ['sub', subject]), 'dir')
        mkdir(fullfile(output_dir, ['sub', subject]));
    end 
    save(fullfile(output_dir, ['sub', subject], ['sub', subject, '_wordRatingContext_run1.mat']), 'onsets_run1', 'durations_run1', 'names');
    save(fullfile(output_dir, ['sub', subject], ['sub', subject, '_wordRatingContext_run2.mat']), 'onsets_run2', 'durations_run2', 'names');
    if add_button
        %% GET TIMING OF BUTTON PRESS
        psychopy_csv_file_run1.button_rts(end) = psychopy_csv_file_run1.button_rts_final(end);
        psychopy_csv_file_run2.button_rts(end) = psychopy_csv_file_run2.button_rts_final(end);
        button_onsets_run1 = cellfun(@(s) str2double(regexp(s, '[\d.]+', 'match')), psychopy_csv_file_run1.button_rts, 'UniformOutput', false);
        button_onsets_run2 = cellfun(@(s) str2double(regexp(s, '[\d.]+', 'match')), psychopy_csv_file_run2.button_rts, 'UniformOutput', false);
        button_onsets_run1 = [button_onsets_run1{:}] - scan_start_run1;
        button_onsets_run2 = [button_onsets_run2{:}] - scan_start_run2;
        % All 0 for button durations (stick function)
        button_durations_run1 = zeros(1, length(button_onsets_run1));
        button_durations_run2 = zeros(1, length(button_onsets_run2));

        %% CORRECT FORMAT FOR SPM: CELL ARRAYS FOR EACH CONDITION
        names = {'word', 'rating', 'context', 'button'};
        onsets_run1 = {stim_onsets_run1, rating_onsets_run1, context_onsets_run1, button_onsets_run1};
        durations_run1 = {stim_durations_run1, rating_durations_run1, context_durations_run1, button_durations_run1};
        onsets_run2 = {stim_onsets_run2, rating_onsets_run2, context_onsets_run2, button_onsets_run2};
        durations_run2 = {stim_durations_run2, rating_durations_run2, context_durations_run2, button_durations_run2};
        % Save the data in SPM format
        if ~exist(fullfile(output_dir, ['sub', subject]), 'dir')
            mkdir(fullfile(output_dir, ['sub', subject]));
        end
        save(fullfile(output_dir, ['sub', subject], ['sub', subject, '_wordRatingContextButton_run1.mat']), 'onsets_run1', 'durations_run1', 'names');
        save(fullfile(output_dir, ['sub', subject], ['sub', subject, '_wordRatingContextButton_run2.mat']), 'onsets_run2', 'durations_run2', 'names');
    end
    
    disp(['Successfully created timing files for subject ', subject]);
end