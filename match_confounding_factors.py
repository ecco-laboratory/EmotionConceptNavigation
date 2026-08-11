import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import matplotlib.lines as mlines
import matplotlib.patches as patches
from sklearn.metrics import pairwise_distances
from sklearn.manifold import MDS
from scipy.stats import spearmanr
from scipy.spatial import procrustes
from scipy.spatial.distance import cosine
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm
import math
import os
#turn off warnings
import warnings
warnings.filterwarnings('ignore')

if_mds = True

cross_validations = ["xModalityRun", "xModality"]
value_factors = ["start_valence_difference", "start_arousal_difference", "end_valence_difference", "end_arousal_difference", "trajectory_length_difference", "valence_displacement_difference", "arousal_displacement_difference", "same_start", "same_end", "same_id", "same_response", 'pattern_similarity']

data_all = []
data_all = []
for cross_validation in cross_validations:
    if if_mds:
        data = pd.read_csv(f"./outputs/singleTrialBetaAnalysis/noICA/incl_all_subs_trials/onoffGridcontrast_multivariate_mdsAlignedcosine_wMoreControls/seedseedAvg/unsmoothed/Subavg/periodicity6/xModality/phi_voxelComponentAverage/singleTrialBeta/csv/design_tbl_va_cat_id_response.csv")
    else:
        data = pd.read_csv(f"./outputs/singleTrialBetaAnalysis/noICA/incl_all_subs_trials/onoffGridcontrast_multivariate_wMoreControls/unsmoothed/Subavg/periodicity6/{cross_validation}/phi_voxelComponentAverage/singleTrialBeta/csv/design_tbl_va_cat_id_response.csv")
    available_factors = [col for col in value_factors if col in data.columns]
    data["cross_validation"] = cross_validation
    data_all.append(data)

data_all = pd.concat(data_all, ignore_index=True)
data_all = data_all[data_all["region"] == "OFC2016ConstantinescuR5"].copy()


def bootstrap_ci(data, target_col, group_vars, num_iter=10000, conf=0.95, random_state=42):
    rng = np.random.default_rng(random_state)
    alpha = (1 - conf) / 2

    results = []
    grouped = data.groupby(group_vars)
    for group_keys, group_df in grouped:
        values = group_df[target_col].dropna().to_numpy()
        if len(values) == 0:
            continue

        #bootstrap resampling
        boot_means = [np.mean(rng.choice(values, size=len(values), replace=True)) for _ in range(num_iter)]
        ci_lower = np.percentile(boot_means, 100 * alpha)
        ci_upper = np.percentile(boot_means, 100 * (1 - alpha))
        mean_val = np.mean(values)

        boot_std = np.std(boot_means)
        d = mean_val / (boot_std*np.sqrt(len(values)))

        if not isinstance(group_keys, tuple):
            group_keys = (group_keys,)

        results.append((*group_keys, mean_val, ci_lower, ci_upper, d))

    columns = group_vars + ["mean", "ci_lower", "ci_upper", "d"]
    return pd.DataFrame(results, columns=columns)

def sign_flip_permutation(data, target_col, group_vars, num_iter=10000, random_state=200, alternative='one-sided',report_bonferroni=False):
    results = []
    grouped = data.groupby(group_vars)
    
    for group_keys, group_df in grouped:
        values = group_df[target_col].dropna().to_numpy()
        if len(values) == 0:
            continue
        
        # Observed test statistic (mean)
        observed_stat = np.mean(values)
        
        # Generate null distribution using sign-flip permutations
        null_stats = []
        n = len(values)
        np.random.seed(random_state)
        
        for _ in range(num_iter):
            # Randomly flip signs
            signs = np.random.choice([-1, 1], size=n)
            permuted_values = values * signs
            null_stat = np.mean(permuted_values)
            null_stats.append(null_stat)
        
        # Calculate one-sided p-value (proportion of null stats >= observed)
        null_stats = np.array(null_stats)
        if alternative == 'one-sided':
            p_val = (np.sum(null_stats >= observed_stat)+1)/(num_iter+1)
        elif alternative == 'two-sided':
            p_val = (np.sum(np.abs(null_stats) >= np.abs(observed_stat))+1)/(num_iter+1)
        results.append((*group_keys, observed_stat, p_val))
        
    if report_bonferroni:
        # get number of groups
        num_groups = len(results)
        results_df = pd.DataFrame(results, columns=group_vars + ["test_stat", "p_val"])
        results_df['bonferroni_p_val'] = results_df['p_val'] * num_groups
        #cap at 1
        results_df['bonferroni_p_val'] = results_df['bonferroni_p_val'].clip(upper=1)
        return results_df
    else:
        return pd.DataFrame(results, columns=group_vars + ["test_stat", "p_val"])

alignment_col, outcome_col = "alignment", "pattern_similarity"
discrete_match_factors = ["same_start", "same_end", "same_id", "same_response"]
continuous_match_factors = [] #['trajectory_length_difference', 'valence_displacement_difference']#"end_valence_difference"]
continuous_all_factors = ["start_valence_difference", "start_arousal_difference", "end_valence_difference", "end_arousal_difference", "trajectory_length_difference", "valence_displacement_difference", "arousal_displacement_difference"]
group_cols = ["region", "subject", "cross_validation", "test_runs"]
n_repeats = 5000
caliper = 0.25
random_state = 42


def exact_and_continuous_match(group_df, rng, caliper=0.25):
    stratum_matches, use_continuous = [], bool(continuous_match_factors)
    for _, stratum_df in group_df.groupby(discrete_match_factors, sort=False, dropna=False):
        aligned = stratum_df.loc[stratum_df[alignment_col] == 1].copy()
        misaligned = stratum_df.loc[stratum_df[alignment_col] == 0].copy()
        if aligned.empty or misaligned.empty: continue
        if not use_continuous:
            n_match = min(len(aligned), len(misaligned))
            matched_aligned_idx = rng.choice(aligned.index.to_numpy(), size=n_match, replace=False)
            matched_misaligned_idx = rng.choice(misaligned.index.to_numpy(), size=n_match, replace=False)
        else:
            combined = pd.concat([aligned[continuous_match_factors], misaligned[continuous_match_factors]], axis=0)
            if combined.isna().any().any(): continue
            scaler = StandardScaler().fit(combined)
            x_aligned = scaler.transform(aligned[continuous_match_factors])
            x_misaligned = scaler.transform(misaligned[continuous_match_factors])
            aligned_order = rng.permutation(len(aligned))
            available_misaligned = list(range(len(misaligned)))
            matched_aligned_idx, matched_misaligned_idx = [], []
            for ai in aligned_order:
                if not available_misaligned: break
                distances = np.linalg.norm(x_misaligned[available_misaligned] - x_aligned[ai], axis=1)
                nearest_position = int(np.argmin(distances))
                nearest_distance = distances[nearest_position]
                if nearest_distance > caliper: continue
                mi = available_misaligned.pop(nearest_position)
                matched_aligned_idx.append(aligned.index[ai])
                matched_misaligned_idx.append(misaligned.index[mi])
        if len(matched_aligned_idx) > 0:
            stratum_matches.extend([group_df.loc[matched_aligned_idx], group_df.loc[matched_misaligned_idx]])
    return pd.concat(stratum_matches, ignore_index=True) if stratum_matches else group_df.iloc[0:0].copy()

rng = np.random.default_rng(random_state)
contrast_rows, count_rows, saved_matched_parts = [], [], []
target_region = "OFC2016ConstantinescuR5"
required_cols = [alignment_col, outcome_col, *discrete_match_factors, *continuous_all_factors, *group_cols]
numeric_cols = [alignment_col, outcome_col, *discrete_match_factors, *continuous_all_factors]

match_data = data_all.loc[data_all["region"] == target_region, required_cols].copy()
match_data[numeric_cols] = match_data[numeric_cols].apply(pd.to_numeric, errors="coerce")
match_data = match_data.dropna(subset=required_cols).copy()
match_data[alignment_col] = match_data[alignment_col].astype(int)
match_data[discrete_match_factors] = match_data[discrete_match_factors].astype(int)

print(f"Rows in match_data: {len(match_data)}")
print(f"Participants: {match_data['subject'].nunique()}")
print(match_data[alignment_col].value_counts())
grouped_data = list(match_data.groupby(group_cols, sort=False))

print(f"Rows available: {len(match_data):,}")
print(f"Fold-level groups: {len(grouped_data):,}")

contrast_rows, count_rows, factor_rows, saved_matched_parts = [], [], [], []

all_factors = discrete_match_factors + continuous_all_factors + ['pattern_similarity'] #continuous_match_factors

for repeat in tqdm(range(n_repeats)):
    for group_keys, group_df in grouped_data:
        matched = exact_and_continuous_match(group_df, rng, caliper)
        if matched.empty: continue
        counts = matched[alignment_col].value_counts()
        if not {0, 1}.issubset(counts.index): continue
        region, subject, cross_validation, test_runs = group_keys
        keys = {"repeat": repeat, "region": region, "subject": subject, "cross_validation": cross_validation, "test_runs": test_runs}
        means = matched.groupby(alignment_col)[outcome_col].mean()
        contrast_rows.append({**keys, "aligned_mean": means.loc[1], "misaligned_mean": means.loc[0], "pattern_similarity_difference": means.loc[1] - means.loc[0]})
        count_rows.append({**keys, "n_aligned": int(counts.loc[1]), "n_misaligned": int(counts.loc[0]), "n_original": len(group_df), "n_matched": len(matched)})
        for factor in all_factors:
            a = matched.loc[matched[alignment_col]==1, factor]
            m = matched.loc[matched[alignment_col]==0, factor]
            pooled_sd = np.sqrt((a.var(ddof=1)+m.var(ddof=1))/2)
            factor_rows.append({**keys, "factor": factor, "factor_type": "discrete" if factor in discrete_match_factors else "continuous", "aligned_mean": a.mean(), "misaligned_mean": m.mean(), "difference": a.mean()-m.mean(), "smd": (a.mean()-m.mean())/pooled_sd if pooled_sd>0 else 0})
        if repeat == 0: saved_matched_parts.append(matched.assign(**keys))

fold_contrasts, matched_counts, factor_balance_fold = pd.DataFrame(contrast_rows), pd.DataFrame(count_rows), pd.DataFrame(factor_rows)
matched_data = pd.concat(saved_matched_parts, ignore_index=True) if saved_matched_parts else pd.DataFrame()
if fold_contrasts.empty: raise RuntimeError("No matched samples were found.")

# Pattern-similarity contrasts
cv_repeat_contrasts = fold_contrasts.groupby(["region", "subject", "cross_validation", "repeat"], as_index=False)[["aligned_mean", "misaligned_mean", "pattern_similarity_difference"]].mean()
subject_repeat_contrasts = cv_repeat_contrasts.groupby(["region", "subject", "repeat"], as_index=False)[["aligned_mean", "misaligned_mean", "pattern_similarity_difference"]].mean()
subject_matched_contrasts = subject_repeat_contrasts.groupby(["region", "subject"], as_index=False).agg(aligned_mean=("aligned_mean", "mean"), misaligned_mean=("misaligned_mean", "mean"), pattern_similarity_difference=("pattern_similarity_difference", "mean"), matching_repeat_sd=("pattern_similarity_difference", "std"))

# Confound balance
factor_balance_cv = factor_balance_fold.groupby(["region", "subject", "cross_validation", "repeat", "factor", "factor_type"], as_index=False)[["aligned_mean", "misaligned_mean", "difference", "smd"]].mean()
factor_balance_subject = factor_balance_cv.groupby(["region", "subject", "repeat", "factor", "factor_type"], as_index=False)[["aligned_mean", "misaligned_mean", "difference", "smd"]].mean()
factor_balance_repeat = factor_balance_subject.groupby(["region", "repeat", "factor", "factor_type"], as_index=False).agg(
    aligned_mean=("aligned_mean", "mean"),
    misaligned_mean=("misaligned_mean", "mean"),
    mean_difference=("difference", "mean"),
    sd_difference=("difference", "std"),
    mean_smd=("smd", "mean"),
    max_abs_smd=("smd", lambda x: x.abs().max()),
    n_subjects=("subject", "nunique"),
)
factor_balance_repeat["se_difference"] = factor_balance_repeat["sd_difference"] / np.sqrt(factor_balance_repeat["n_subjects"])
factor_balance_across_repeats = factor_balance_repeat.groupby(["region", "factor", "factor_type"], as_index=False).agg(
    mean_aligned=("aligned_mean", "mean"),
    mean_misaligned=("misaligned_mean", "mean"),
    mean_difference=("mean_difference", "mean"),
    sd_difference_across_repeats=("mean_difference", "std"),
    min_difference=("mean_difference", "min"),
    max_difference=("mean_difference", "max"),
    mean_abs_smd=("mean_smd", lambda x: x.abs().mean()),
    max_abs_smd=("max_abs_smd", "max"),
)

# Retention and inference
matched_counts["retained_proportion"] = matched_counts["n_matched"] / matched_counts["n_original"]
print("\nMatching retention:")
print(matched_counts["retained_proportion"].describe()[["mean", "min", "max"]].round(3))
print("\nBalance across all matching repetitions:")
print(factor_balance_across_repeats.round(6).to_string(index=False))

matched_ci = bootstrap_ci(subject_matched_contrasts, "pattern_similarity_difference", ["region"], 10000, 0.95, 42)
matched_test = sign_flip_permutation(subject_matched_contrasts, "pattern_similarity_difference", ["region"], 10000, 42, "two-sided", False)
matched_stats = matched_ci.merge(matched_test, on="region")

# Save
if if_mds:
    output_dir = "./outputs/matched_confounding_factors_analysis_mds/5000repeats"
else:
    output_dir = "./outputs/matched_confounding_factors_analysis/5000repeats"
os.makedirs(output_dir, exist_ok=True)
prefix = "OFC2016ConstantinescuR5_6fold_subavg_exact_continuous_matching"

subject_repeat_contrasts.to_csv(f"{output_dir}/{prefix}_subject_contrasts_by_repeat.csv", index=False)
fold_matching_results = fold_contrasts.merge(matched_counts, on=["repeat", "region", "subject", "cross_validation", "test_runs"], how="left", validate="one_to_one")
fold_matching_results.to_csv(f"{output_dir}/{prefix}_fold_results_by_repeat.csv", index=False)
factor_balance_subject.to_csv(f"{output_dir}/{prefix}_confound_balance_subject_by_repeat.csv", index=False)
factor_balance_repeat.to_csv(f"{output_dir}/{prefix}_confound_balance_group_by_repeat.csv", index=False)
factor_balance_across_repeats.to_csv(f"{output_dir}/{prefix}_confound_balance_across_repeats.csv", index=False)
matched_data.to_csv(f"{output_dir}/{prefix}_balance_repeat_0.csv", index=False)
subject_matched_contrasts.to_csv(f"{output_dir}/{prefix}_subject_contrasts_averaged.csv", index=False)
matched_stats.to_csv(f"{output_dir}/{prefix}_group_statistics.csv", index=False)