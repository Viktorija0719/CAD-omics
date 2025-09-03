#!/usr/bin/env python3
"""
DESCRIPTION:
    Statistical comparison of clinical variables between two groups (status 0 vs 1).

INPUTS:
    (1) Clinical data CSV file
        - Rows = patients, Columns = clinical features (continuous + categorical).
        - Must contain a 'sample' column for patient IDs.
    (2) Status CSV file
        - Two columns: {sample, status}, where status ∈ {0,1}.
    (3) Output CSV path
        - Destination file where the statistical results will be written.

ANALYSIS:
    - Continuous variables:
        * Shapiro–Wilk test is used to check normality in each group.
        * If both groups are normally distributed → Welch's t-test is performed.
          Summary statistics are reported as mean ± SD.
        * If at least one group is non-normal → Mann–Whitney U test is performed.
          Summary statistics are reported as median (IQR).
    - Categorical variables:
        * 2×2 contingency tables → Fisher’s exact test.
        * Larger tables → Chi-squared test.
        * Summary shows group-wise counts and percentages.

OUTPUT:
    A CSV file with columns:
        - Variable : clinical variable name
        - Test     : statistical test applied
        - p-value  : computed p-value
        - Summary  : descriptive stats per group (continuous) 
                     or contingency counts (categorical)

USAGE:
    python stats_compare.py \
        --clinical data/processed/clinical/stability_meta_missing.csv \
        --status   data/processed/status/stability_status.csv \
        --output   results/stability_statistical_comparison.csv

DEPENDENCIES:
    pandas, numpy, scipy
"""

import argparse
import pandas as pd
import numpy as np
from scipy.stats import shapiro, ttest_ind, mannwhitneyu, fisher_exact, chi2_contingency

continuous_columns = [
    "PRS", "Age", "Weight", "Height", "BMI", "Pack_years", "WBC", "RBC", "HGB", "HCT", "PLT",
    "NE_abs", "LY_abs", "MO_abs", "Eo_Abs", "Total_cholesterol", "LDL", "HDL", "TG",
    "Creatinine", "ALAT", "ASAT", "Glucose", "HbA1c", "LAD_stenosis",
    "LCX_stenosis", "RCA_stenosis", "Target_fibrotic", "Target_lipidic", "Target_necrotic",
    "Target_Calcific", "Plaque_necrolipidic_tissue", "Dist_fibrotic", "Dist_lipidic", "Dist_necrotic",
    "Dist_calcific"
]

categorical_columns = [
    "Sex", "Smoking_1", "Smoking_2", "Positive_family_history",
    "Art_hipert", "Congestive_heart_failure",
    "Previous_PCI", "Previous_MI"
]

def fmt_summary(series: pd.Series, is_normal: bool) -> str:
    x = series.dropna().to_numpy()
    if x.size == 0:
        return "NA"
    if is_normal:
        return f"{np.mean(x):.2f} \u00B1 {np.std(x, ddof=1):.2f}"
    q1, q3 = np.percentile(x, [25, 75])
    iqr = q3 - q1
    med = np.median(x)
    return f"{med:.2f} (IQR: {iqr:.2f})"

def compare_continuous(df, col):
    g0 = df.loc[df["status"] == 0, col].dropna()
    g1 = df.loc[df["status"] == 1, col].dropna()

    if len(g0) < 3 or len(g1) < 3:
        return {"Variable": col, "Test": "NA", "p-value": np.nan,
                "Summary": "Insufficient data"}

    # Normality per group (Shapiro; skip for very large n)
    _, p0 = shapiro(g0) if len(g0) < 5000 else (None, 1.0)
    _, p1 = shapiro(g1) if len(g1) < 5000 else (None, 1.0)
    normal0 = p0 > 0.05
    normal1 = p1 > 0.05

    if normal0 and normal1:
        stat, pval = ttest_ind(g0, g1, equal_var=False)
        test = "t-test"
        sum0 = fmt_summary(g0, True)
        sum1 = fmt_summary(g1, True)
    else:
        stat, pval = mannwhitneyu(g0, g1, alternative="two-sided")
        test = "Mann-Whitney U"
        sum0 = fmt_summary(g0, False)
        sum1 = fmt_summary(g1, False)

    summary = f"status 0: {sum0} | status 1: {sum1}"
    return {"Variable": col, "Test": test, "p-value": float(pval), "Summary": summary}

def compare_categorical(df, col):
    # Ensure categorical/binary-like columns are treated as categorical
    contingency = pd.crosstab(df["status"], df[col])
    if contingency.empty or contingency.shape[0] != 2:
        return {"Variable": col, "Test": "NA", "p-value": np.nan, "Summary": "Insufficient data"}

    if contingency.shape == (2, 2):
        _, pval = fisher_exact(contingency)
        test = "Fisher's exact"
    else:
        try:
            _, pval, _, _ = chi2_contingency(contingency)
            test = "Chi-sq"
        except Exception:
            return {"Variable": col, "Test": "NA", "p-value": np.nan, "Summary": "Insufficient data"}

    # Compact summary of counts per level
    # e.g., {"0": {"No": 10, "Yes": 5}, "1": {"No": 7, "Yes": 12}}
    summary = contingency.to_dict()
    return {"Variable": col, "Test": test, "p-value": float(pval), "Summary": summary}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--clinical", required=True, help="Path to clinical CSV")
    parser.add_argument("--status", required=True, help="Path to status CSV")
    parser.add_argument("--output", required=True, help="Path to output CSV with results")
    args = parser.parse_args()

    clin = pd.read_csv(args.clinical)
    status = pd.read_csv(args.status)

    # Ensure 'status' is numeric 0/1
    status["status"] = pd.to_numeric(status["status"], errors="coerce")

    df = clin.merge(status, on="sample", how="inner")

    results = []

    # Continuous
    for col in continuous_columns:
        if col in df.columns:
            results.append(compare_continuous(df, col))

    # Categorical
    for col in categorical_columns:
        if col in df.columns:
            results.append(compare_categorical(df, col))

    res_df = pd.DataFrame(results, columns=["Variable", "Test", "p-value", "Summary"])
    res_df.to_csv(args.output, index=False)
    print(f"✅ Results saved to {args.output}")

if __name__ == "__main__":
    main()
