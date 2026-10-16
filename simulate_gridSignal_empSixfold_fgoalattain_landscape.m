
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

emp_rA = 0.0801;
emp_rM = 0.0729;
emp_raw_mean = 0.0765;
emp_delta = 0.0072;
emp_sd_raw_mid = 0.0262;
emp_sd_delta = 0.009;
emp_sd_rA = 0.0262;
emp_corr_rArM = 0.936;

n_grid = 10;
n_mc = 500;                
rng_seed_land = 42;         

sim_mode = getenv('SIM_MODE');
if isempty(sim_mode), sim_mode = 'run'; end

cache_base = fullfile(project_dir, 'outputs', 'simulations', 'gridSignal_empSixfold');
if normalize_B, cache_base = [cache_base, '_normalize_B']; end
if normalize_g, cache_base = [cache_base, '_normalize_g']; end
cache_file = fullfile(cache_base, 'empirical_empSixfold_cache.mat');

fit_dir = fullfile(project_dir, 'outputs', 'simulations', 'gridSignal_empSixfold_fgoalattain');
if normalize_B, fit_dir = [fit_dir, '_normalize_B']; end
if normalize_g, fit_dir = [fit_dir, '_normalize_g']; end
opt_result_path = fullfile(fit_dir, 'fgoalattain_opt_result.mat');

out_dir = fullfile(fit_dir, 'landscape');
shard_dir = fullfile(out_dir, 'shards');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
if ~exist(shard_dir, 'dir'), mkdir(shard_dir); end

fprintf('Landscape | normalize_B=%d normalize_g=%d | n_grid=%d n_mc=%d\n', normalize_B, normalize_g, n_grid, n_mc);
fprintf('Fit dir: %s\n', fit_dir);
fprintf('Out dir: %s\n', out_dir);

assert(exist(opt_result_path, 'file') == 2, 'Missing optimum: %s\nRun fgoalattain first.', opt_result_path);
Sopt = load(opt_result_path);
b0_hat = Sopt.b0_hat;
sigma_b_hat = Sopt.sigma_b_hat;
A_hat = Sopt.A_hat;
fprintf('x_hat: b0=%g  sigma_b=%g  A=%g\n', b0_hat, sigma_b_hat, A_hat);

% Local grids: n_grid points centered on hats with fixed step sizes
b0_step = 0.05;
sigma_b_step = 0.02;
A_step = 0.04;
offs = ((0:n_grid-1) - (n_grid-1)/2);   % e.g. n_grid=10 -> -4.5:1:4.5

b0_grid = max(1e-4, b0_hat + offs * b0_step);
sigma_b_grid = max(0, sigma_b_hat + offs * sigma_b_step);
A_grid = max(0, A_hat + offs * A_step);

% Nearest grid points to hats (for optional plot lines)
b0_hat_snap = b0_grid(argmin_abs(b0_grid, b0_hat));
sigma_b_hat_snap = sigma_b_grid(argmin_abs(sigma_b_grid, sigma_b_hat));
A_hat_snap = A_grid(argmin_abs(A_grid, A_hat));
fprintf('Grids (hat-centered, fixed step; n_grid=%d):\n', n_grid);
fprintf('  b0:      step=%g  [%g, %g]  (hat=%g, snap=%g)\n', b0_step, b0_grid(1), b0_grid(end), b0_hat, b0_hat_snap);
fprintf('  sigma_b: step=%g  [%g, %g]  (hat=%g, snap=%g)\n', sigma_b_step, sigma_b_grid(1), sigma_b_grid(end), sigma_b_hat, sigma_b_hat_snap);
fprintf('  A:       step=%g  [%g, %g]  (hat=%g, snap=%g)\n', A_step, A_grid(1), A_grid(end), A_hat, A_hat_snap);

if strcmpi(sim_mode, 'aggregate')
    aggregate_landscape_shards(out_dir, shard_dir, b0_hat, sigma_b_hat, A_hat, ...
        b0_grid, sigma_b_grid, A_grid, emp_rA, emp_rM, emp_raw_mean, emp_delta, ...
        emp_sd_raw_mid, emp_sd_delta, emp_corr_rArM);
    return;
end

%% Load empirical structure
fprintf('Loading angles + empirical sixfold cache...\n');
[angle_data_all, subject_ids, modality_run_ids] = load_and_prepare_angle_data(subjects, runs, spm1stlevel_dir, theta_dir, theta_source);
n_subjects = numel(subjects);
n_trials_per_sub = zeros(n_subjects, 1);
for s = 1:n_subjects
    n_trials_per_sub(s) = sum(strcmp(subject_ids, subjects{s}));
end
[B_emp, g_emp, phi_emp_rad, L_resid, ~] = load_empirical_empSixfold( ...
    subjects, spm1stlevel_dir, roi_mask, n_voxels, normalize_B, normalize_g, ...
    ridge_cov, periodicity, angle_data_all, subject_ids, cache_file);

S_cache = load(cache_file);
sd_eps = nan(n_subjects, n_voxels);
for s = 1:n_subjects
    sd_eps(s, :) = sqrt(diag(S_cache.Sigma_resid_cell{s}))';
end
inv_sd_mean = mean(1 ./ sd_eps(:), 'omitnan');

%% Common random numbers (identical for every task / grid point)
fprintf('Drawing eps for landscape (%d reps, seed_base=%d)...\n', n_mc, rng_seed_land);
z_b_all = cell(n_mc, 1);
z_eps_all = cell(n_mc, 1);
for rep = 1:n_mc
    rng(rng_seed_land + 100000 + rep);
    z_b_all{rep} = randn(n_subjects, 1);
    z_eps_all{rep} = cell(n_subjects, 1);
    for s = 1:n_subjects
        z_eps_all{rep}{s} = randn(n_trials_per_sub(s), n_voxels);
    end
end

%% Run
task_id = str2double(getenv('SLURM_ARRAY_TASK_ID'));
n_tasks = 3 * n_grid;
if isnan(task_id) || task_id < 1 || task_id > n_tasks
    fprintf('No/invalid SLURM_ARRAY_TASK_ID; running ALL %d tasks serially.\n', n_tasks);
    task_list = 1:n_tasks;
else
    task_list = task_id;
    fprintf('SLURM task %d / %d\n', task_id, n_tasks);
end


meta_path = fullfile(out_dir, 'landscape_meta.mat');
if ~exist(meta_path, 'file') || ~isnan(task_id) && task_id == 1
    save(meta_path, 'b0_hat', 'sigma_b_hat', 'A_hat', ...
        'b0_hat_snap', 'sigma_b_hat_snap', 'A_hat_snap', ...
        'b0_grid', 'sigma_b_grid', 'A_grid', 'b0_step', 'sigma_b_step', 'A_step', ...
        'n_grid', 'n_mc', 'rng_seed_land', 'normalize_B', 'normalize_g', 'inv_sd_mean');
    meta_tbl = table(b0_hat, sigma_b_hat, A_hat, b0_hat_snap, sigma_b_hat_snap, A_hat_snap, ...
        n_grid, n_mc, rng_seed_land, normalize_B, normalize_g, ...
        'VariableNames', {'b0_hat', 'sigma_b_hat', 'A_hat', ...
        'b0_hat_snap', 'sigma_b_hat_snap', 'A_hat_snap', ...
        'n_grid', 'n_mc', 'rng_seed_land', 'normalize_B', 'normalize_g'});
    writetable(meta_tbl, fullfile(out_dir, 'landscape_meta.csv'));
end

for tid = task_list(:)'
    plane = ceil(tid / n_grid);
    idx = mod(tid - 1, n_grid) + 1;
    run_landscape_task(plane, idx, b0_hat, sigma_b_hat, A_hat, ...
        b0_grid, sigma_b_grid, A_grid, n_mc, z_b_all, z_eps_all, L_resid, ...
        B_emp, g_emp, phi_emp_rad, angle_data_all, subject_ids, modality_run_ids, ...
        subjects, periodicity, cross_validations, inv_sd_mean, shard_dir, ...
        emp_rA, emp_rM, emp_raw_mean, emp_delta, emp_sd_raw_mid, emp_sd_delta, emp_corr_rArM);
end

fprintf('Done task(s).\n');


%% ===================== landscape task / aggregate ============================

function run_landscape_task(plane, idx, b0_hat, sigma_b_hat, A_hat, ...
    b0_grid, sigma_b_grid, A_grid, n_mc, z_b_all, z_eps_all, L_resid, ...
    B_emp, g_emp, phi_emp_rad, angle_data_all, subject_ids, modality_run_ids, ...
    subjects, periodicity, cross_validations, inv_sd_mean, shard_dir, ...
    emp_rA, emp_rM, emp_raw_mean, emp_delta, emp_sd_raw_mid, emp_sd_delta, emp_corr)

    n_subjects = numel(subjects);
    switch plane
        case 1  % fix b0, vary sigma_b (outer) x A (inner)
            b0 = b0_hat;
            sigma_b = sigma_b_grid(idx);
            inner_vals = A_grid;
            inner_name = 'A';
            tag = sprintf('plane1_fixb0_sigmab%02d', idx);
            fprintf('Plane 1 | b0=hat=%g | sigma_b[%d]=%g | sweep A (%d)\n', ...
                b0, idx, sigma_b, numel(inner_vals));
        case 2  % fix sigma_b, vary b0 (outer) x A (inner)
            b0 = b0_grid(idx);
            sigma_b = sigma_b_hat;
            inner_vals = A_grid;
            inner_name = 'A';
            tag = sprintf('plane2_fixsb_b0%02d', idx);
            fprintf('Plane 2 | sigma_b=hat=%g | b0[%d]=%g | sweep A (%d)\n', ...
                sigma_b, idx, b0, numel(inner_vals));
        case 3  % fix A, vary b0 (outer) x sigma_b (inner)
            b0 = b0_grid(idx);
            A_fixed = A_hat;
            inner_vals = sigma_b_grid;
            inner_name = 'sigma_b';
            tag = sprintf('plane3_fixA_b0%02d', idx);
            fprintf('Plane 3 | A=hat=%g | b0[%d]=%g | sweep sigma_b (%d)\n', ...
                A_fixed, idx, b0, numel(inner_vals));
        otherwise
            error('Unknown plane %d', plane);
    end

    rep_rows = {};
    for rep = 1:n_mc
        if mod(rep, 50) == 1 || rep == n_mc
            fprintf('  %s rep %d/%d\n', tag, rep, n_mc);
        end
        z_b = z_b_all{rep};
        z_eps_cell = z_eps_all{rep};

        for j = 1:numel(inner_vals)
            if plane == 1
                A = inner_vals(j);
                b_s = b0 + sigma_b * z_b;
            elseif plane == 2
                A = inner_vals(j);
                b_s = b0 + sigma_b * z_b;
            else
                sigma_b = inner_vals(j);
                A = A_fixed;
                b_s = b0 + sigma_b * z_b;
            end

            eps_sim_cell = cell(n_subjects, 1);
            for s = 1:n_subjects
                eps_sim_cell{s} = z_eps_cell{s} * L_resid{s}';
            end

            [sub_tbl, ~] = run_full_analysis(A, B_emp, g_emp, b_s, ...
                phi_emp_rad, eps_sim_cell, angle_data_all, subject_ids, modality_run_ids, ...
                subjects, periodicity, cross_validations);
            pooled = pool_cvs_within_subject(sub_tbl);
            raw_mid = 0.5 * (pooled.r_aligned + pooled.r_misaligned);
            SNR = A * inv_sd_mean;

            rep_rows{end+1} = {rep, plane, b0, sigma_b, A, ...
                mean(pooled.r_aligned, 'omitnan'), mean(pooled.r_misaligned, 'omitnan'), ...
                mean(raw_mid, 'omitnan'), mean(pooled.delta_r, 'omitnan'), ...
                std(pooled.r_aligned, 'omitnan'), std(pooled.r_misaligned, 'omitnan'), ...
                std(raw_mid, 'omitnan'), std(pooled.delta_r, 'omitnan'), ...
                corr(pooled.r_aligned, pooled.r_misaligned, 'rows', 'complete'), SNR};
        end
    end

    rep_tbl = cell2table(vertcat(rep_rows{:}), 'VariableNames', { ...
        'rep', 'plane', 'b0', 'sigma_b', 'A', ...
        'mean_rA', 'mean_rM', 'mean_raw_mid', 'mean_delta', ...
        'sd_subject_rA', 'sd_subject_rM', 'sd_subject_raw_mid', 'sd_subject_delta', ...
        'corr_rA_rM', 'SNR'});
    writetable(rep_tbl, fullfile(shard_dir, sprintf('rep_metrics_%s.csv', tag)));

    summary_tbl = summarize_landscape_reps(rep_tbl, emp_rA, emp_rM, emp_raw_mean, ...
        emp_delta, emp_sd_raw_mid, emp_sd_delta, emp_corr);
    writetable(summary_tbl, fullfile(shard_dir, sprintf('mc_summary_%s.csv', tag)));
    fprintf('Wrote shards for %s (%d rows summary)\n', tag, height(summary_tbl));
end

function summary_tbl = summarize_landscape_reps(rep_tbl, emp_rA, emp_rM, emp_raw_mean, ...
    emp_delta, emp_sd_raw_mid, emp_sd_delta, emp_corr)
    keys = unique(rep_tbl(:, {'plane', 'b0', 'sigma_b', 'A'}), 'rows');
    metrics = {'mean_rA', 'mean_rM', 'mean_raw_mid', 'mean_delta', ...
        'sd_subject_rA', 'sd_subject_rM', 'sd_subject_raw_mid', 'sd_subject_delta', ...
        'corr_rA_rM', 'SNR'};
    rows = {};
    for i = 1:height(keys)
        mask = abs(rep_tbl.plane - keys.plane(i)) < 1e-12 & ...
            abs(rep_tbl.b0 - keys.b0(i)) < 1e-12 & ...
            abs(rep_tbl.sigma_b - keys.sigma_b(i)) < 1e-12 & ...
            abs(rep_tbl.A - keys.A(i)) < 1e-12;
        block = rep_tbl(mask, :);
        for m = 1:numel(metrics)
            name = metrics{m};
            x = block.(name);
            rows{end+1} = {keys.plane(i), keys.b0(i), keys.sigma_b(i), keys.A(i), name, ...
                mean(x, 'omitnan'), median(x, 'omitnan'), ...
                prctile(x, 2.5), prctile(x, 97.5), std(x, 'omitnan'), ...
                emp_rA, emp_rM, emp_raw_mean, emp_delta, emp_sd_raw_mid, emp_sd_delta, emp_corr};
        end
    end
    summary_tbl = cell2table(vertcat(rows{:}), 'VariableNames', { ...
        'plane', 'b0', 'sigma_b', 'A', 'metric', 'mc_mean', 'mc_median', ...
        'si95_lo', 'si95_hi', 'mc_sd', ...
        'emp_rA', 'emp_rM', 'emp_raw_mean', 'emp_delta', 'emp_sd_raw_mid', 'emp_sd_delta', 'emp_corr_rArM'});
end

function aggregate_landscape_shards(out_dir, shard_dir, b0_hat, sigma_b_hat, A_hat, ...
    b0_grid, sigma_b_grid, A_grid, emp_rA, emp_rM, emp_raw_mean, emp_delta, ...
    emp_sd_raw_mid, emp_sd_delta, emp_corr)
    fprintf('Aggregating landscape shards in %s ...\n', shard_dir);
    sum_files = dir(fullfile(shard_dir, 'mc_summary_plane*.csv'));
    assert(~isempty(sum_files), 'No mc_summary_plane*.csv in %s', shard_dir);
    cells = cell(numel(sum_files), 1);
    for i = 1:numel(sum_files)
        cells{i} = readtable(fullfile(shard_dir, sum_files(i).name));
    end
    summary = vertcat(cells{:});
    writetable(summary, fullfile(out_dir, 'landscape_mc_summary.csv'));

    rep_files = dir(fullfile(shard_dir, 'rep_metrics_plane*.csv'));
    if ~isempty(rep_files)
        cells = cell(numel(rep_files), 1);
        for i = 1:numel(rep_files)
            cells{i} = readtable(fullfile(shard_dir, rep_files(i).name));
        end
        writetable(vertcat(cells{:}), fullfile(out_dir, 'landscape_rep_metrics_all.csv'));
    end

    % Long-form convenience tables per plane (for notebook line plots)
    for plane = 1:3
        sub = summary(summary.plane == plane, :);
        writetable(sub, fullfile(out_dir, sprintf('landscape_mc_summary_plane%d.csv', plane)));
    end

    fprintf('Wrote landscape_mc_summary.csv (%d rows) and per-plane CSVs\n', height(summary));
    fprintf('Hats: b0=%g sigma_b=%g A=%g\n', b0_hat, sigma_b_hat, A_hat);
end

function ix = argmin_abs(grid, val)
    [~, ix] = min(abs(grid - val));
end


%% ===================== pipeline helpers ================================

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
