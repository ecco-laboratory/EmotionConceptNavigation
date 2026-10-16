clear; clc;

[project_dir, ~, ~, theta_dir, ~, ~, ~, ~, ~, ~, ~, spm1stlevel_dir, ~, ~, subjects] = set_up_dirs_constants();
spm1stlevel_dir = fullfile([spm1stlevel_dir, '_noICA_wButton'], 'unsmoothed');

n_voxels = 69;
periodicity = 6;
theta_source = 'Subavg';
runs = [1, 2];
cross_validations = {'xModalityRun', 'xModality'};

roi_mask = fullfile(project_dir, 'masks', 'OFC_2016Constantinescu_r5.nii');
normalize_B = true;
normalize_g = true;
ridge_cov = 1e-8;

% Empirical targets
emp_rA = 0.0801;
emp_rM = 0.0729;
emp_raw_mean = 0.0765;          % (rA + rM) / 2
emp_delta = 0.0072;
emp_sd_raw_mid = 0.0262;        
emp_sd_delta = 0.009;
emp_sd_rA = 0.0262;
emp_corr_rArM = 0.936;

% Optimization / MC settings
n_opt_reps = 40;                % reps inside each objective eval
n_mc = 500;                     % independent post-fit MC
rng_seed_opt = 42;
rng_seed_mc = 900000;           % different from optimization seeds
x0 = [0.5; 0.1; 0.3];        % [b0; sigma_b; A]
lb = [1e-4; 0; 0];
ub = [2.0; 2.0; 2.0];

sim_mode = getenv('SIM_MODE');
if isempty(sim_mode), sim_mode = 'run'; end

cache_base = fullfile(project_dir, 'outputs', 'simulations', 'gridSignal_empSixfold');
if normalize_B, cache_base = [cache_base, '_normalize_B']; end
if normalize_g, cache_base = [cache_base, '_normalize_g']; end
cache_file = fullfile(cache_base, 'empirical_empSixfold_cache.mat');
if ~exist(cache_base, 'dir'), mkdir(cache_base); end

out_dir = fullfile(project_dir, 'outputs', 'simulations', 'gridSignal_empSixfold_fgoalattain');
if normalize_B, out_dir = [out_dir, '_normalize_B']; end
if normalize_g, out_dir = [out_dir, '_normalize_g']; end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
fprintf('fgoalattain empSixfold | normalize_B=%d | normalize_g=%d\n', normalize_B, normalize_g);
fprintf('Cache: %s\n', cache_file);
fprintf('Output: %s\n', out_dir);

%% Load angles + empirical structure (cache if available)
fprintf('Loading empirical %s angles for %d subjects...\n', theta_source, numel(subjects));
[angle_data_all, subject_ids, modality_run_ids] = load_and_prepare_angle_data(subjects, runs, spm1stlevel_dir, theta_dir, theta_source);
n_subjects = numel(subjects);
n_trials_per_sub = zeros(n_subjects, 1);
for s = 1:n_subjects
    n_trials_per_sub(s) = sum(strcmp(subject_ids, subjects{s}));
end

fprintf('Loading / fitting empirical sixfold (B, g, phi, Sigma_resid)...\n');
[B_emp, g_emp, phi_emp_rad, L_resid, diag_tbl] = load_empirical_empSixfold(subjects, spm1stlevel_dir, roi_mask, n_voxels, normalize_B, normalize_g, ridge_cov, periodicity, angle_data_all, subject_ids, cache_file);
writetable(diag_tbl, fullfile(out_dir, 'empSixfold_diagnostics.csv'));

% SNR = mean_{s,v}(A / sd_eps(s,v)) = A * mean(1./sd_eps)
assert(exist(cache_file, 'file') == 2, 'Expected cache with Sigma_resid after load_empirical');
S_cache = load(cache_file);
assert(isfield(S_cache, 'Sigma_resid_cell'), 'Cache missing Sigma_resid_cell: %s', cache_file);
sd_eps = nan(n_subjects, n_voxels);
for s = 1:n_subjects
    sd_eps(s, :) = sqrt(diag(S_cache.Sigma_resid_cell{s}))';
end
inv_sd_mean = mean(1 ./ sd_eps(:), 'omitnan');

%% Shared draws of epsilon for optimization (fixed across all objective evals)
fprintf('Drawing common random numbers for optimization (%d reps)...\n', n_opt_reps);
z_b_opt = cell(n_opt_reps, 1);
z_eps_opt = cell(n_opt_reps, 1);
for rep = 1:n_opt_reps
    rng(rng_seed_opt + 100000 + rep);
    z_b_opt{rep} = randn(n_subjects, 1);
    z_eps_opt{rep} = cell(n_subjects, 1);
    for s = 1:n_subjects
        z_eps_opt{rep}{s} = randn(n_trials_per_sub(s), n_voxels);
    end
end

opt_hist_path = fullfile(out_dir, 'fgoalattain_eval_history.csv');
opt_result_path = fullfile(out_dir, 'fgoalattain_opt_result.mat');

do_opt = strcmpi(sim_mode, 'run') || strcmpi(sim_mode, 'opt_only');
do_mc = strcmpi(sim_mode, 'run') || strcmpi(sim_mode, 'mc_only');

if do_opt
    % Fresh history file
    if exist(opt_hist_path, 'file'), delete(opt_hist_path); end
    hist_hdr = {'eval_id', 'b0', 'sigma_b', 'A', ...
        'mean_rA', 'mean_rM', 'mean_raw_mid', 'mean_delta', ...
        'sd_subject_rA', 'sd_subject_rM', 'sd_subject_raw_mid', 'sd_subject_delta', ...
        'corr_rA_rM', 'SNR', ...
        'e_raw_mean', 'e_delta', 'e_raw_SD', 'e_delta_SD', 'e_sum'};
    writetable(cell2table(cell(0, numel(hist_hdr)), 'VariableNames', hist_hdr), opt_hist_path);

    eval_count = 0;

    emp_vec = [emp_raw_mean; emp_delta; emp_sd_raw_mid; emp_sd_delta];
    goal = zeros(4, 1);
    weight = [1 2 1 1];

    objfun = @(x) goalattain_objective(x, emp_vec, inv_sd_mean, ...
        n_opt_reps, z_b_opt, z_eps_opt, L_resid, ...
        B_emp, g_emp, phi_emp_rad, angle_data_all, subject_ids, modality_run_ids, ...
        subjects, periodicity, cross_validations, opt_hist_path);

    % Wrapper to count evals 
    opts = optimoptions('fgoalattain', ...
        'Display', 'iter', ...
        'MaxFunctionEvaluations', 200, ...
        'MaxIterations', 80, ...
        'FunctionTolerance', 1e-6, ...
        'StepTolerance', 1e-4, ...
        'UseParallel', false);

    fprintf('\n=== fgoalattain: x=[b0, sigma_b, A], goals=0, equal weights ===\n');
    fprintf('x0 = [%g, %g, %g]\n', x0(1), x0(2), x0(3));
    [x_hat, fval, attainfactor, exitflag, output] = fgoalattain(objfun, x0, goal, weight, [], [], [], [], lb, ub, [], opts);

    b0_hat = x_hat(1); sigma_b_hat = x_hat(2); A_hat = x_hat(3);
    fprintf('\nBest fit: b0=%g  sigma_b=%g  A=%g\n', b0_hat, sigma_b_hat, A_hat);
    fprintf('F (squared abs err) = [%g, %g, %g, %g]  attainfactor=%g  exitflag=%d\n', fval(1), fval(2), fval(3), fval(4), attainfactor, exitflag);

    % Final stats at optimum
    [stats_hat, e_hat] = simulate_summary_crn(x_hat, n_opt_reps, z_b_opt, z_eps_opt, L_resid, B_emp, g_emp, phi_emp_rad, angle_data_all, subject_ids, modality_run_ids, subjects, periodicity, cross_validations, inv_sd_mean, emp_vec);
    emp_names = {'raw_mean', 'delta', 'raw_SD', 'delta_SD'};

    save(opt_result_path, 'x_hat', 'b0_hat', 'sigma_b_hat', 'A_hat', 'fval', 'attainfactor', ...
        'exitflag', 'output', 'x0', 'lb', 'ub', 'goal', 'weight', 'emp_vec', ...
        'emp_rA', 'emp_rM', 'emp_raw_mean', 'emp_delta', 'emp_sd_raw_mid', 'emp_sd_delta', ...
        'stats_hat', 'e_hat', 'n_opt_reps', 'rng_seed_opt', 'normalize_B', 'normalize_g', ...
        'inv_sd_mean');

    best_tbl = table(b0_hat, sigma_b_hat, A_hat, ...
        stats_hat.mean_rA, stats_hat.mean_rM, stats_hat.mean_raw_mid, stats_hat.mean_delta, ...
        stats_hat.sd_subject_rA, stats_hat.sd_subject_rM, stats_hat.sd_subject_raw_mid, stats_hat.sd_subject_delta, ...
        stats_hat.corr_rA_rM, stats_hat.SNR, ...
        e_hat(1), e_hat(2), e_hat(3), e_hat(4), sum(e_hat), attainfactor, exitflag, ...
        'VariableNames', {'b0', 'sigma_b', 'A', ...
        'mean_rA', 'mean_rM', 'mean_raw_mid', 'mean_delta', ...
        'sd_subject_rA', 'sd_subject_rM', 'sd_subject_raw_mid', 'sd_subject_delta', ...
        'corr_rA_rM', 'SNR', ...
        'e_raw_mean', 'e_delta', 'e_raw_SD', 'e_delta_SD', 'e_sum', 'attainfactor', 'exitflag'});
    writetable(best_tbl, fullfile(out_dir, 'fgoalattain_best_fit.csv'));
    fprintf('Wrote %s\n', opt_result_path);
else
    assert(exist(opt_result_path, 'file') == 2, 'Missing %s (run opt first)', opt_result_path);
    Sopt = load(opt_result_path);
    b0_hat = Sopt.b0_hat; sigma_b_hat = Sopt.sigma_b_hat; A_hat = Sopt.A_hat;
    inv_sd_mean = Sopt.inv_sd_mean;
    fprintf('Loaded prior optimum: b0=%g sigma_b=%g A=%g\n', b0_hat, sigma_b_hat, A_hat);
end

%% post-fit replicates
if do_mc
    fprintf('\n=== Post-fit replicates: n_mc=%d, seed_base=%d ===\n', n_mc, rng_seed_mc);
    rep_rows = cell(n_mc, 1);
    subj_rows = {};
    for rep = 1:n_mc
        if mod(rep, 25) == 1 || rep == n_mc
            fprintf('  replicate %d / %d\n', rep, n_mc);
        end
        rng(rng_seed_mc + rep);
        z_b = randn(n_subjects, 1);
        b_s = b0_hat + sigma_b_hat * z_b;
        eps_sim_cell = cell(n_subjects, 1);
        for s = 1:n_subjects
            z = randn(n_trials_per_sub(s), n_voxels);
            eps_sim_cell{s} = z * L_resid{s}';
        end
        [sub_tbl, ~] = run_full_analysis(A_hat, B_emp, g_emp, b_s, ...
            phi_emp_rad, eps_sim_cell, angle_data_all, subject_ids, modality_run_ids, ...
            subjects, periodicity, cross_validations);
        pooled = pool_cvs_within_subject(sub_tbl);
        raw_mid = 0.5 * (pooled.r_aligned + pooled.r_misaligned);
        SNR = A_hat * inv_sd_mean;
        rep_rows{rep} = {rep, b0_hat, sigma_b_hat, A_hat, ...
            mean(pooled.r_aligned, 'omitnan'), mean(pooled.r_misaligned, 'omitnan'), ...
            mean(raw_mid, 'omitnan'), mean(pooled.delta_r, 'omitnan'), ...
            std(pooled.r_aligned, 'omitnan'), std(pooled.r_misaligned, 'omitnan'), ...
            std(raw_mid, 'omitnan'), std(pooled.delta_r, 'omitnan'), ...
            corr(pooled.r_aligned, pooled.r_misaligned, 'rows', 'complete'), SNR};
        for s = 1:height(pooled)
            subj_rows{end+1} = {rep, b0_hat, sigma_b_hat, A_hat, pooled.subject{s}, ...
                pooled.b_s(s), pooled.r_aligned(s), pooled.r_misaligned(s), ...
                0.5 * (pooled.r_aligned(s) + pooled.r_misaligned(s)), pooled.delta_r(s), ...
                pooled.phi_err_deg(s)};
        end
    end

    rep_tbl = cell2table(vertcat(rep_rows{:}), 'VariableNames', { ...
        'rep', 'b0', 'sigma_b', 'A', ...
        'mean_rA', 'mean_rM', 'mean_raw_mid', 'mean_delta', ...
        'sd_subject_rA', 'sd_subject_rM', 'sd_subject_raw_mid', 'sd_subject_delta', ...
        'corr_rA_rM', 'SNR'});
    subj_tbl = cell2table(vertcat(subj_rows{:}), 'VariableNames', { ...
        'rep', 'b0', 'sigma_b', 'A', 'subject', 'b_s', ...
        'rA', 'rM', 'raw_mid', 'delta', 'phi_err_deg'});

    writetable(rep_tbl, fullfile(out_dir, 'postfit_rep_metrics.csv'));
    writetable(subj_tbl, fullfile(out_dir, 'postfit_subject_level.csv'));

    metric_names = {'mean_rA', 'mean_rM', 'mean_raw_mid', 'mean_delta', ...
        'sd_subject_rA', 'sd_subject_rM', 'sd_subject_raw_mid', 'sd_subject_delta', ...
        'corr_rA_rM', 'SNR'};
    sum_rows = cell(numel(metric_names), 1);
    for m = 1:numel(metric_names)
        name = metric_names{m};
        x = rep_tbl.(name);
        sum_rows{m} = {b0_hat, sigma_b_hat, A_hat, name, ...
            mean(x, 'omitnan'), median(x, 'omitnan'), ...
            prctile(x, 2.5), prctile(x, 97.5), std(x, 'omitnan'), ...
            emp_rA, emp_rM, emp_raw_mean, emp_delta, emp_sd_raw_mid, emp_sd_delta, emp_corr_rArM};
    end
    summary_tbl = cell2table(vertcat(sum_rows{:}), 'VariableNames', { ...
        'b0', 'sigma_b', 'A', 'metric', 'mc_mean', 'mc_median', 'si95_lo', 'si95_hi', 'mc_sd', ...
        'emp_rA', 'emp_rM', 'emp_raw_mean', 'emp_delta', 'emp_sd_raw_mid', 'emp_sd_delta', 'emp_corr_rArM'});
    writetable(summary_tbl, fullfile(out_dir, 'postfit_mc_summary.csv'));
    fprintf('Wrote post-fit MC summary: %s\n', fullfile(out_dir, 'postfit_mc_summary.csv'));
    disp(summary_tbl(:, {'metric', 'mc_mean', 'si95_lo', 'si95_hi'}));
end

fprintf('Done.\n');



%% ===================== optimization helpers ==================================

function F = goalattain_objective(x, emp_vec, inv_sd_mean, ...
    n_opt_reps, z_b_opt, z_eps_opt, L_resid, ...
    B_emp, g_emp, phi_emp_rad, angle_data_all, subject_ids, modality_run_ids, ...
    subjects, periodicity, cross_validations, hist_path)

    x = x(:);
    [stats, e] = simulate_summary_crn(x, n_opt_reps, z_b_opt, z_eps_opt, L_resid, ...
        B_emp, g_emp, phi_emp_rad, angle_data_all, subject_ids, modality_run_ids, ...
        subjects, periodicity, cross_validations, inv_sd_mean, emp_vec);
    F = e(:);

    row = table(NaN, x(1), x(2), x(3), ...
        stats.mean_rA, stats.mean_rM, stats.mean_raw_mid, stats.mean_delta, ...
        stats.sd_subject_rA, stats.sd_subject_rM, stats.sd_subject_raw_mid, stats.sd_subject_delta, ...
        stats.corr_rA_rM, stats.SNR, ...
        e(1), e(2), e(3), e(4), sum(e), ...
        'VariableNames', {'eval_id', 'b0', 'sigma_b', 'A', ...
        'mean_rA', 'mean_rM', 'mean_raw_mid', 'mean_delta', ...
        'sd_subject_rA', 'sd_subject_rM', 'sd_subject_raw_mid', 'sd_subject_delta', ...
        'corr_rA_rM', 'SNR', ...
        'e_raw_mean', 'e_delta', 'e_raw_SD', 'e_delta_SD', 'e_sum'});
    if exist(hist_path, 'file')
        try
            old = readtable(hist_path);
            if height(old) == 0
                row.eval_id = 1;
                writetable(row, hist_path);
            else
                row.eval_id = height(old) + 1;
                writetable([old; row], hist_path);
            end
        catch
            row.eval_id = 1;
            writetable(row, hist_path);
        end
    else
        row.eval_id = 1;
        writetable(row, hist_path);
    end
    fprintf('  eval %d | b0=%.4f sb=%.4f A=%.4f | raw=%.4f d=%.4f sd_raw=%.4f sd_d=%.4f | e=[%.3g %.3g %.3g %.3g]\n', ...
        row.eval_id, x(1), x(2), x(3), stats.mean_raw_mid, stats.mean_delta, ...
        stats.sd_subject_raw_mid, stats.sd_subject_delta, e(1), e(2), e(3), e(4));
end

function [stats, e] = simulate_summary_crn(x, n_reps, z_b_opt, z_eps_opt, L_resid, ...
    B_emp, g_emp, phi_emp_rad, angle_data_all, subject_ids, modality_run_ids, ...
    subjects, periodicity, cross_validations, inv_sd_mean, emp_vec)

    b0 = max(x(1), 1e-8);
    sigma_b = max(x(2), 0);
    A = max(x(3), 0);
    n_subjects = numel(subjects);

    mean_rA = nan(n_reps, 1);
    mean_rM = nan(n_reps, 1);
    mean_raw = nan(n_reps, 1);
    mean_delta = nan(n_reps, 1);
    sd_rA = nan(n_reps, 1);
    sd_rM = nan(n_reps, 1);
    sd_raw = nan(n_reps, 1);
    sd_delta = nan(n_reps, 1);
    corr_AM = nan(n_reps, 1);

    for rep = 1:n_reps
        b_s = b0 + sigma_b * z_b_opt{rep};
        eps_sim_cell = cell(n_subjects, 1);
        for s = 1:n_subjects
            eps_sim_cell{s} = z_eps_opt{rep}{s} * L_resid{s}';
        end
        [sub_tbl, ~] = run_full_analysis(A, B_emp, g_emp, b_s, ...
            phi_emp_rad, eps_sim_cell, angle_data_all, subject_ids, modality_run_ids, ...
            subjects, periodicity, cross_validations);
        pooled = pool_cvs_within_subject(sub_tbl);
        raw_mid = 0.5 * (pooled.r_aligned + pooled.r_misaligned);
        mean_rA(rep) = mean(pooled.r_aligned, 'omitnan');
        mean_rM(rep) = mean(pooled.r_misaligned, 'omitnan');
        mean_raw(rep) = mean(raw_mid, 'omitnan');
        mean_delta(rep) = mean(pooled.delta_r, 'omitnan');
        sd_rA(rep) = std(pooled.r_aligned, 'omitnan');
        sd_rM(rep) = std(pooled.r_misaligned, 'omitnan');
        sd_raw(rep) = std(raw_mid, 'omitnan');
        sd_delta(rep) = std(pooled.delta_r, 'omitnan');
        corr_AM(rep) = corr(pooled.r_aligned, pooled.r_misaligned, 'rows', 'complete');
    end

    stats = struct();
    stats.mean_rA = mean(mean_rA, 'omitnan');
    stats.mean_rM = mean(mean_rM, 'omitnan');
    stats.mean_raw_mid = mean(mean_raw, 'omitnan');
    stats.mean_delta = mean(mean_delta, 'omitnan');
    stats.sd_subject_rA = mean(sd_rA, 'omitnan');
    stats.sd_subject_rM = mean(sd_rM, 'omitnan');
    stats.sd_subject_raw_mid = mean(sd_raw, 'omitnan');
    stats.sd_subject_delta = mean(sd_delta, 'omitnan');
    stats.corr_rA_rM = mean(corr_AM, 'omitnan');
    stats.SNR = A * inv_sd_mean;

    sim_vec = [stats.mean_raw_mid; stats.mean_delta; stats.sd_subject_raw_mid; stats.sd_subject_delta];
    e = (sim_vec - emp_vec).^2;  % squared error
end


%% ===================== other funcs ================

function [B_emp, g_emp, phi_emp_rad, L_resid, diag_tbl] = load_empirical_empSixfold(subjects, spm1stlevel_dir, roi_mask, n_voxels, normalize_B, normalize_g, ridge_cov, periodicity, angle_data_all, subject_ids, cache_file)

    cache_ok = false;
    if exist(cache_file, 'file')
        try
            S = load(cache_file);
            if isfield(S, 'normalize_B') && S.normalize_B == normalize_B && ...
                    isfield(S, 'normalize_g') && S.normalize_g == normalize_g && ...
                    isfield(S, 'roi_mask') && strcmp(S.roi_mask, roi_mask) && ...
                    isfield(S, 'B_emp') && size(S.B_emp, 2) == n_voxels && ...
                    isfield(S, 'subjects') && isequal(S.subjects(:), subjects(:)) && ...
                    isfield(S, 'g_emp') && isfield(S, 'phi_emp_rad') && isfield(S, 'L_resid')
                B_emp = S.B_emp;
                g_emp = S.g_emp;
                phi_emp_rad = S.phi_emp_rad;
                L_resid = S.L_resid;
                if isfield(S, 'diag_tbl')
                    diag_tbl = S.diag_tbl;
                else
                    diag_tbl = table();
                end
                cache_ok = true;
                fprintf('  loaded empSixfold cache: %s\n', cache_file);
            end
        catch ME
            warning('%s', ['Ignoring bad cache; recomputing. ', ME.message]);
            try, delete(cache_file); end
        end
    end
    if cache_ok
        if isempty(diag_tbl) || height(diag_tbl) == 0
        end
        return;
    end

    n_subjects = numel(subjects);
    B_emp = nan(n_subjects, n_voxels);
    g_emp = nan(n_subjects, n_voxels);
    phi_emp_rad = nan(n_subjects, 1);
    L_resid = cell(n_subjects, 1);
    Sigma_resid_cell = cell(n_subjects, 1);
    beta_cos_emp_all = nan(n_subjects, n_voxels);
    beta_sin_emp_all = nan(n_subjects, n_voxels);

    diag_rows = cell(n_subjects, 1);

    for s = 1:n_subjects
        sub = subjects{s};
        angles = angle_data_all(strcmp(subject_ids, sub));
        angles = angles(:);
        n_trials = numel(angles);

        singleTrial_dir = fullfile(spm1stlevel_dir, ['sub', sub], 'singleTrial');
        load(fullfile(singleTrial_dir, 'SPM.mat'));
        stim_idx = find(contains(SPM.xX.name, 'trial'));
        if numel(stim_idx) ~= n_trials
            error('Subject %s: n_stim_betas=%d but n_angles=%d', sub, numel(stim_idx), n_trials);
        end
        beta_files = cell(numel(stim_idx), 1);
        for i = 1:numel(stim_idx)
            beta_files{i} = fullfile(singleTrial_dir, sprintf('beta_%04d.nii', stim_idx(i)));
        end

        beta_obj = apply_mask(fmri_data(beta_files), roi_mask);
        beta_mat = double(beta_obj.dat)';  % trials x voxels
        if size(beta_mat, 2) ~= n_voxels
            error('Subject %s: expected %d voxels, got %d', sub, n_voxels, size(beta_mat, 2));
        end
        if size(beta_mat, 1) ~= n_trials
            error('Subject %s: beta trials %d != angle trials %d', sub, size(beta_mat, 1), n_trials);
        end

        % Stable pattern B_s
        B = mean(beta_mat, 1)';
        B = B - mean(B);
        if normalize_B
            sdB = std(B);
            if sdB < eps, error('Subject %s: B has zero SD', sub); end
            B = B / sdB;
        end
        B_emp(s, :) = B';

        % Voxelwise quadrature: beta = alpha + c*cos(6*theta) + s*sin(6*theta)
        X_noPhi = [cos(periodicity * deg2rad(angles)), sin(periodicity * deg2rad(angles))];
        X_design = [ones(n_trials, 1), X_noPhi];
        beta_cos_emp = nan(n_voxels, 1);
        beta_sin_emp = nan(n_voxels, 1);
        alpha_emp = nan(n_voxels, 1);
        E_sixfold_resid = nan(n_trials, n_voxels);
        for v = 1:n_voxels
            bb = X_design \ beta_mat(:, v);
            alpha_emp(v) = bb(1);
            beta_cos_emp(v) = bb(2);
            beta_sin_emp(v) = bb(3);
            yhat = X_design * bb;
            E_sixfold_resid(:, v) = beta_mat(:, v) - yhat;
        end
        beta_cos_emp_all(s, :) = beta_cos_emp';
        beta_sin_emp_all(s, :) = beta_sin_emp';

        % ROI-level empirical orientation (voxel-component average)
        phi_emp = atan2(mean(beta_sin_emp), mean(beta_cos_emp)) / periodicity;
        phi_emp_rad(s) = phi_emp;

        % Voxelwise sixfold amplitudes
        g = sqrt(beta_cos_emp.^2 + beta_sin_emp.^2);
        g_mean_raw = mean(g);
        g_rms_raw = sqrt(mean(g.^2));
        if normalize_g
            if g_mean_raw < eps
                error('Subject %s: mean(g_emp) too small to normalize', sub);
            end
            g = g / g_mean_raw;  % mean(g)=1 after scaling
        end
        g_emp(s, :) = g';

        % Variance before/after sixfold removal (trial-demeaned for fair compare)
        E_demean = beta_mat - mean(beta_mat, 1);
        var_before = mean(E_demean(:).^2);
        var_after = mean(E_sixfold_resid(:).^2);
        frac_var_removed = 1 - var_after / max(var_before, eps);

        % Sigma_resid from sixfold-residualized trial patterns
        % use cholesky decomposition to make episilon_st to have empirical covariance matrix
        E = E_sixfold_resid - mean(E_sixfold_resid, 1);
        Sigma_resid = cov(E, 1);
        Sigma_resid = (Sigma_resid + Sigma_resid') / 2;
        Sigma_resid = Sigma_resid + ridge_cov * eye(n_voxels);
        [L, pflag] = chol(Sigma_resid, 'lower');
        if pflag > 0
            ridge = ridge_cov;
            for k = 1:8
                ridge = ridge * 10;
                Sigma_try = Sigma_resid - ridge_cov * eye(n_voxels) + ridge * eye(n_voxels);
                [L, pflag] = chol(Sigma_try, 'lower');
                if pflag == 0
                    Sigma_resid = Sigma_try;
                    break;
                end
            end
            if pflag > 0
                error('Subject %s: failed to Cholesky-factor Sigma_resid', sub);
            end
        end
        Sigma_resid_cell{s} = Sigma_resid;
        L_resid{s} = L;

        % corr of upper-tri cov entries demeaned vs residualized
        Sigma_demean = cov(E_demean, 1);
        Sigma_demean = (Sigma_demean + Sigma_demean') / 2;
        ut = triu(true(n_voxels), 1);
        corr_cov = corr(Sigma_demean(ut), Sigma_resid(ut), 'rows', 'complete');

        R_s = sqrt(mean(beta_cos_emp)^2 + mean(beta_sin_emp)^2);

        diag_rows{s} = {sub, rad2deg(phi_emp), mean(g), std(g), min(g), max(g), ...
            g_mean_raw, g_rms_raw, R_s, var_before, var_after, frac_var_removed, ...
            trace(Sigma_resid), corr_cov, normalize_g};

        fprintf(['  sub%s: phi=%f deg | g mean/sd/min/max=%g/%g/%g/%g', ...
            ' | R_s=%g | frac_var_removed=%f | trace(Sigma_resid)=%g\n'], ...
            sub, rad2deg(phi_emp), mean(g), std(g), min(g), max(g), R_s, ...
            frac_var_removed, trace(Sigma_resid));
    end

    diag_tbl = cell2table(vertcat(diag_rows{:}), 'VariableNames', { ...
        'subject', 'phi_emp_deg', 'g_mean', 'g_sd', 'g_min', 'g_max', ...
        'g_mean_raw', 'g_rms_raw', 'R_s', 'var_before_sixfold', 'var_after_sixfold', ...
        'frac_var_removed', 'trace_Sigma_resid', 'corr_cov_demean_vs_resid', 'normalize_g'});

    try
        save(cache_file, 'B_emp', 'g_emp', 'phi_emp_rad', 'L_resid', 'Sigma_resid_cell', ...
            'beta_cos_emp_all', 'beta_sin_emp_all', 'diag_tbl', ...
            'normalize_B', 'normalize_g', 'roi_mask', 'subjects', '-v7.3');
        fprintf('  saved empSixfold cache: %s\n', cache_file);
    catch ME
        warning('%s', ['Cache save skipped (likely array race): ', ME.message]);
    end
end


function xyz_mm = get_roi_xyz_mm(masked_obj)
    % After apply_mask, dat has n_roi voxels; volInfo.xyzlist stays brain-sized.
    % Non-removed voxels match dat row order.
    vi = masked_obj.volInfo;
    n = size(masked_obj.dat, 1);
    if isempty(masked_obj.removed_voxels) || numel(masked_obj.removed_voxels) ~= size(vi.xyzlist, 1)
        error('get_roi_xyz_mm: unexpected removed_voxels / xyzlist sizes');
    end
    kept = ~masked_obj.removed_voxels;
    if sum(kept) ~= n
        error('get_roi_xyz_mm: kept=%d but dat has %d voxels', sum(kept), n);
    end
    xyz_vox = double(vi.xyzlist(kept, :));
    xyz_mm = (vi.mat * [xyz_vox(:, 1:3), ones(size(xyz_vox, 1), 1)]')';
    xyz_mm = xyz_mm(:, 1:3);
end

function [sub_tbl, fold_tbl] = run_full_analysis(A, B_all, g_all, b_s, ...
    phi_true_rad, eps_cell, angle_data_all, subject_ids, modality_run_ids, subjects, ...
    periodicity, cross_validations)
    fold_rows = {};
    sub_rows = {};
    for s = 1:numel(subjects)
        sub = subjects{s};
        mask = strcmp(subject_ids, sub);
        angles = angle_data_all(mask);
        mod_ids = modality_run_ids(mask);
        n_trials = numel(angles);
        B = B_all(s, :)';
        g = g_all(s, :)';
        phi_true = phi_true_rad(s);
        b = b_s(s);
        eps_s = eps_cell{s};
        if size(eps_s, 1) ~= n_trials
            error('eps size mismatch for subject %s (eps %d, trials %d)', sub, size(eps_s, 1), n_trials);
        end
        cos_term = cos(periodicity * deg2rad(angles - rad2deg(phi_true)));
        % beta = b*B + A*g*cos[6(theta-phi)] + eps_sim
        betas = (b * B) * ones(1, n_trials) + A * (g * cos_term') + eps_s';
        betas = betas';

        for cv = 1:numel(cross_validations)
            cv_name = cross_validations{cv};
            [train_folds, test_folds] = get_cv_folds(cv_name);
            rA_folds = nan(numel(train_folds), 1);
            rM_folds = nan(numel(train_folds), 1);
            d_folds = nan(numel(train_folds), 1);
            phi_err_folds = nan(numel(train_folds), 1);
            nA_folds = nan(numel(train_folds), 1);
            nM_folds = nan(numel(train_folds), 1);
            for f = 1:numel(train_folds)
                train_mask = ismember(mod_ids, train_folds{f});
                test_mask = ismember(mod_ids, test_folds{f});
                phi_hat = estimate_phi_voxel_component_average(betas(train_mask, :), angles(train_mask), periodicity);
                phi_err_folds(f) = circular_abs_error_deg(phi_hat, phi_true, periodicity);
                [delta, rA, rM, nA, nM] = perform_multivariate_contrast_with_means(betas(test_mask, :), angles(test_mask), phi_hat, periodicity);
                rA_folds(f) = rA; rM_folds(f) = rM; d_folds(f) = delta;
                nA_folds(f) = nA; nM_folds(f) = nM;
                fold_rows{end+1} = {sub, cv_name, A, b, rad2deg(phi_true), rad2deg(phi_hat), phi_err_folds(f), rA, rM, delta, nA, nM};
            end
            sub_rows{end+1} = {sub, cv_name, A, b, rad2deg(phi_true), ...
                mean(rA_folds, 'omitnan'), mean(rM_folds, 'omitnan'), mean(d_folds, 'omitnan'), ...
                mean(phi_err_folds, 'omitnan'), mean(nA_folds, 'omitnan'), mean(nM_folds, 'omitnan')};
        end
    end
    sub_tbl = cell2table(vertcat(sub_rows{:}), 'VariableNames', ...
        {'subject', 'cv', 'A', 'b_s', 'phi_true_deg', ...
         'r_aligned', 'r_misaligned', 'delta_r', 'phi_err_deg', 'n_aligned_pairs', 'n_misaligned_pairs'});
    if ~isempty(fold_rows)
        fold_tbl = cell2table(vertcat(fold_rows{:}), 'VariableNames', ...
            {'subject', 'cv', 'A', 'b_s', 'phi_true_deg', 'phi_hat_deg', ...
             'phi_err_deg', 'r_aligned', 'r_misaligned', 'delta_r', 'n_aligned_pairs', 'n_misaligned_pairs'});
    else
        fold_tbl = table();
    end
end

function [train_folds, test_folds] = get_cv_folds(cv_name)
    if strcmp(cv_name, 'xModalityRun')
        train_folds = {{'face1', 'word1'}, {'face1', 'word2'}, {'face2', 'word1'}, {'face2', 'word2'}};
        test_folds  = {{'face2', 'word2'}, {'face2', 'word1'}, {'face1', 'word2'}, {'face1', 'word1'}};
    elseif strcmp(cv_name, 'xModality')
        train_folds = {{'face1', 'face2'}, {'word1', 'word2'}};
        test_folds  = {{'word1', 'word2'}, {'face1', 'face2'}};
    else
        error('Unknown CV: %s', cv_name);
    end
end

function phi = estimate_phi_voxel_component_average(Y_train, X_train, periodicity)
    X_train = X_train(:);
    X_noPhi = [cos(periodicity * deg2rad(X_train)), sin(periodicity * deg2rad(X_train))];
    X_design = [ones(size(X_noPhi, 1), 1), X_noPhi];
    n_vox = size(Y_train, 2);
    beta_cos = nan(n_vox, 1); beta_sin = nan(n_vox, 1);
    for v = 1:n_vox
        bb = X_design \ Y_train(:, v);
        beta_cos(v) = bb(2); beta_sin(v) = bb(3);
    end
    phi = atan2(mean(beta_sin), mean(beta_cos)) / periodicity;
end

function err_deg = circular_abs_error_deg(phi_hat, phi_true, periodicity)
    period_deg = 360 / periodicity;
    d = abs(rad2deg(phi_hat) - rad2deg(phi_true));
    d = mod(d, period_deg);
    err_deg = min(d, period_deg - d);
end

function [contrast_value, aligned_mean, misaligned_mean, n_aligned, n_misaligned] = perform_multivariate_contrast_with_means(singleTrial_betas, angle_data, phi_value, current_periodicity)
    phi_deg = rad2deg(phi_value);
    period_deg = 360 / current_periodicity;
    bin_half = (period_deg / 2) / 2;
    align_half = bin_half / 2;
    angle_data = angle_data(:);
    nTrials = size(singleTrial_betas, 1);
    if nTrials < 2
        contrast_value = NaN; aligned_mean = NaN; misaligned_mean = NaN;
        n_aligned = 0; n_misaligned = 0; return;
    end
    angle_modP = mod(angle_data - phi_deg, period_deg);
    aligned_idx = find((abs(angle_modP) < align_half) | (abs(angle_modP - period_deg) < align_half));
    if isempty(aligned_idx)
        contrast_value = NaN; aligned_mean = NaN; misaligned_mean = NaN;
        n_aligned = 0; n_misaligned = 0; return;
    end
    R = corr(singleTrial_betas');
    R(1:nTrials+1:end) = NaN;
    aligned_corrs = []; misaligned_corrs = [];
    for i = 1:numel(aligned_idx)
        t = aligned_idx(i);
        d_on_period = mod(abs(angle_data(t) - angle_data), period_deg);
        diff_to_0 = min(d_on_period, period_deg - d_on_period);
        diff_to_half = abs(d_on_period - period_deg / 2);
        diff_to_half = min(diff_to_half, period_deg - diff_to_half);
        is_aligned_pair = diff_to_0 < bin_half; is_aligned_pair(t) = false;
        is_misaligned_pair = diff_to_half < bin_half; is_misaligned_pair(t) = false;
        aligned_corrs = [aligned_corrs; R(t, is_aligned_pair)'];
        misaligned_corrs = [misaligned_corrs; R(t, is_misaligned_pair)'];
    end
    n_aligned = numel(aligned_corrs); n_misaligned = numel(misaligned_corrs);
    if n_aligned == 0 || n_misaligned == 0
        aligned_mean = NaN; misaligned_mean = NaN; contrast_value = NaN;
    else
        aligned_mean = mean(aligned_corrs, 'omitnan');
        misaligned_mean = mean(misaligned_corrs, 'omitnan');
        contrast_value = aligned_mean - misaligned_mean;
    end
end

function pooled = pool_cvs_within_subject(sub_tbl)
    subs = unique(sub_tbl.subject, 'stable');
    rows = cell(numel(subs), 1);
    for i = 1:numel(subs)
        idx = strcmp(sub_tbl.subject, subs{i});
        rows{i} = {subs{i}, mean(sub_tbl.A(idx)), mean(sub_tbl.b_s(idx)), ...
            mean(sub_tbl.phi_true_deg(idx)), mean(sub_tbl.r_aligned(idx), 'omitnan'), ...
            mean(sub_tbl.r_misaligned(idx), 'omitnan'), mean(sub_tbl.delta_r(idx), 'omitnan'), ...
            mean(sub_tbl.phi_err_deg(idx), 'omitnan'), mean(sub_tbl.n_aligned_pairs(idx), 'omitnan'), ...
            mean(sub_tbl.n_misaligned_pairs(idx), 'omitnan')};
    end
    pooled = cell2table(vertcat(rows{:}), 'VariableNames', ...
        {'subject', 'A', 'b_s', 'phi_true_deg', 'r_aligned', 'r_misaligned', ...
         'delta_r', 'phi_err_deg', 'n_aligned_pairs', 'n_misaligned_pairs'});
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
        first_stim_trial = spm_regressor_names(find(stim_trials, 1, 'first'));
        last_stim_trial = spm_regressor_names(find(stim_trials, 1, 'last'));
        if contains(first_stim_trial, 'face') && contains(last_stim_trial, 'word')
            modalities = {'face', 'word'};
        elseif contains(first_stim_trial, 'word') && contains(last_stim_trial, 'face')
            modalities = {'word', 'face'};
        else
            error('Modality order unclear for subject %s.', current_subject);
        end
        angles_all_modalities_runs = [];
        angles_modality_run_idx = {};
        for m = 1:length(modalities)
            current_modality = modalities{m};
            for r = 1:length(runs)
                angle_file = fullfile(theta_dir, ['sub', current_subject], ['sub', current_subject, '_', current_modality, '_thetas_run', num2str(runs(r)), '.mat']);
                load(angle_file);
                angle_data = thetas.(current_angle_source)';
                angle_avg_data = thetas.Subavg';
                if strcmp(current_angle_source, 'Subspec')
                    nan_idx = isnan(angle_data);
                    angle_data(nan_idx) = angle_avg_data(nan_idx);
                end
                angles_all_modalities_runs = [angles_all_modalities_runs; angle_data];
                current_modality_run_idx = [current_modality, num2str(runs(r))];
                angles_modality_run_idx = [angles_modality_run_idx; ...
                    repmat({current_modality_run_idx}, size(angle_data, 1), 1)];
            end
        end
        angle_data_all = [angle_data_all; angles_all_modalities_runs];
        subject_ids = [subject_ids; repmat({current_subject}, size(angles_all_modalities_runs, 1), 1)];
        modality_run_ids = [modality_run_ids; angles_modality_run_idx];
    end
end
