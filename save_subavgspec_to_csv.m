[project_dir, fmri_data_dir, ~, theta_dir, category_dir, ~, ~, ~, spm_timing_dir, motion_regressors_dir, nuisance_regressors_dir, spm1stlevel_dir, ~, phi_dir, subjects] = set_up_dirs_constants();

modalities = {'face', 'word'};
runs = [1, 2];
thetas_avg_allsubs_runs = [];
thetas_spec_allsubs_runs = [];
all_subjects = [];
all_modalities = [];
all_runs = [];


for s = 1:length(subjects)
    sub_id = subjects{s};
    for m = 1:length(modalities)
        modality = modalities{m};
        for r = 1:length(runs)
            run = runs(r);

            %get thetas from subavg and subspec
            theta_file = fullfile(theta_dir, ['sub', sub_id], ['sub', sub_id, '_', modality, '_thetas_run', num2str(run), '.mat']);
            load(theta_file);
            thetas_subavg = thetas.Subavg';
            thetas_subspec = thetas.Subspec';
            thetas_avg_allsubs_runs = [thetas_avg_allsubs_runs; thetas_subavg];
            thetas_spec_allsubs_runs = [thetas_spec_allsubs_runs; thetas_subspec];
            all_subjects = [all_subjects; repmat(string(sub_id), size(thetas_subavg, 1), 1)];
            all_modalities = [all_modalities; repmat(string(modality), size(thetas_subavg, 1), 1)];
            all_runs = [all_runs; repmat(string(run), size(thetas_subavg, 1), 1)];
        end
    end
end

thetas_avg_allsubs_runs = array2table([thetas_avg_allsubs_runs, thetas_spec_allsubs_runs, all_subjects, all_modalities, all_runs], ...
                                'VariableNames', {'theta_avg', 'theta_spec', 'subject', 'modality', 'run'});
writetable(thetas_avg_allsubs_runs, fullfile(project_dir, 'outputs', 'thetas_subavgspec.csv'));

