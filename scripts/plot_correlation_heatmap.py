#!/usr/bin/env python3
"""
DESCRIPTION:
    Plot a correlation heatmap (miRNA vs clinical variables) from the CSV produced by the
    miRNA–clinical correlation script. Supports optional filtering by method, effect size,
    and FDR, plus friendly renaming of clinical variable labels and miRNA IDs.

INPUT:
    --input   Path to correlation CSV with columns:
              ['Variable1','Variable2','Correlation Coefficient','p.value','p.adjusted','Method']
              where Variable1 = clinical variable, Variable2 = miRNA ID.
    --output  Path to save the heatmap image (e.g., PNG/PDF/SVG)

OPTIONS:
    --method        Filter by correlation method: 'both' (default), 'kendall', 'point-biserial'
    --abs-threshold Keep |Correlation Coefficient| > this value (default: 0.20)
    --fdr           Keep rows with p.adjusted < this value (default: 0.05)

USAGE:
    docker run -it --rm \
  -u $(id -u):$(id -g) \
  -v "$(pwd)":/app \
  cad-omics \
  python scripts/plot_correlation_heatmap.py \
      --input results/correlation/stability_miRNA_corr_combined_results.csv \
      --output results/correlation/stability_miRNA_corr_heatmap.png \
      --method both --abs-threshold 0.20 --fdr 0.05

DEPENDENCIES:
    pandas, matplotlib, seaborn
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap


def parse_args():
    p = argparse.ArgumentParser(description="Plot correlation heatmap from correlation CSV results.")
    p.add_argument("--input", required=True, help="Path to the input correlation CSV file")
    p.add_argument("--output", required=True, help="Path to save the output heatmap image")
    p.add_argument("--method", choices=["both", "kendall", "point-biserial"], default="both",
                   help="Correlation method to include [default: both]")
    p.add_argument("--abs-threshold", type=float, default=0.20,
                   help="Absolute correlation threshold [default: 0.20]")
    p.add_argument("--fdr", type=float, default=0.05,
                   help="Adjusted p-value (FDR) cutoff [default: 0.05]")
    return p.parse_args()


def main():
    args = parse_args()

    # Load correlation results
    df = pd.read_csv(args.input)

    # Standardize method filter
    method_map = {
        "kendall": "Kendall",
        "point-biserial": "Point-Biserial",
        "both": None,
    }
    wanted_method = method_map[args.method]
    if wanted_method is not None and "Method" in df.columns:
        df = df[df["Method"].str.casefold() == wanted_method.casefold()].copy()

    # Rename clinical variables for readability
    rename_map = {
        'Age': 'Age', 'Sex': 'Sex', 'BMI': 'BMI', 'WBC': 'WBC', 'RBC': 'RBC', 'HGB': 'HGB', 'HCT': 'HCT', 'PLT': 'PLT',
        'NE_abs': 'NE abs', 'LY_abs': 'LY abs', 'MO_abs': 'MO abs', 'Eo_Abs': 'Eo Abs',
        'Total_cholesterol': 'Total Cholesterol',
        'LDL': 'LDL', 'HDL': 'HDL', 'TG': 'TG', 'Creatinine': 'Creatinine', 'ALAT': 'ALAT', 'ASAT': 'ASAT',
        'Glucose': 'Glucose', 'HbA1c': 'HbA1c',
        'Target_fibrotic': 'Fibrotic tissue', 'Target_lipidic': 'Lipidic tissue',
        'Target_necrotic': 'Necrotic tissue', 'Target_Calcific': 'Calcific tissue',
        'Smoking_1': 'Smoking', 'Smoking_2': 'Quitted smoking',
        'Positive_family_history': 'Positive family history',
        'Art_hipert': 'Arterial Hypertension',
        'Previous_PCI': 'Previous PCI', 'Previous_MI': 'Previous MI'
    }
    if "Variable1" in df.columns:
        df["Variable1"] = df["Variable1"].replace(rename_map)

    # Clean miRNA labels: strip 'hsa-' prefix to show 'miR-...'
    def clean_mirna(x: str) -> str:
        if isinstance(x, str) and x.startswith("hsa-"):
            return x.replace("hsa-", "", 1)
        return x

    if "Variable2" in df.columns:
        df["Variable2"] = df["Variable2"].apply(clean_mirna)

    # # Ensure required columns exist
    # needed = {"Variable1", "Variable2", "Correlation Coefficient", "p.adjusted"}
    # missing = needed - set(df.columns)
    # if missing:
    #     raise ValueError(f"Missing required column(s) in input CSV: {', '.join(sorted(missing))}")

    # # Filter by |r| and FDR
    # df = df[np.isfinite(df["Correlation Coefficient"])].copy()
    # filtered = df[(df["p.adjusted"] < args.fdr) &
    #               (df["Correlation Coefficient"].abs() > args.abs_threshold)].copy()
    

    # Ensure required columns exist
    needed = {"Variable1", "Variable2", "Correlation Coefficient", "p.value"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"Missing required column(s) in input CSV: {', '.join(sorted(missing))}")

    # Filter by |r| and FDR
    df = df[np.isfinite(df["Correlation Coefficient"])].copy()
    filtered = df[(df["p.value"] < args.fdr) &
                  (df["Correlation Coefficient"].abs() > args.abs_threshold)].copy()


    if filtered.empty:
        print("No correlations passed the filters; saving an empty heatmap placeholder.")
        # Create a minimal empty figure rather than erroring out
        plt.figure(figsize=(8, 3))
        plt.text(0.5, 0.5, "No correlations passed filters", ha="center", va="center")
        plt.axis("off")
        plt.savefig(args.output, dpi=300, bbox_inches="tight")
        return

    # Pivot for heatmap: rows = miRNAs, columns = clinical variables
    heatmap_data = filtered.pivot(index="Variable2", columns="Variable1", values="Correlation Coefficient")

    # Retain only desired clinical variable columns in a sensible order (based on rename_map)
    desired_order = [v for v in rename_map.values() if v in heatmap_data.columns]
    # Keep any extra columns (e.g., PRS) at the end
    extras = [c for c in heatmap_data.columns if c not in desired_order]
    ordered_cols = desired_order + sorted(extras)
    heatmap_data = heatmap_data.reindex(columns=ordered_cols)

    # Custom diverging colormap (blue → white → red)
    custom_cmap = LinearSegmentedColormap.from_list("custom_colormap", ["#8e82fe", "white", "#ef4026"])

    # Plot
    plt.figure(figsize=(max(10, len(ordered_cols) * 0.5), max(8, len(heatmap_data) * 0.3)))
    sns.heatmap(
        heatmap_data,
        cmap=custom_cmap,
        vmin=-1, vmax=1,
        annot=False,
        linewidths=0.5,
        cbar_kws={'label': 'Correlation Coefficient (r)'},
        square=False
    )
    title_bits = [f"Method: {wanted_method or 'Both'}", f"|r|>{args.abs_threshold}", f"FDR<{args.fdr}"]
    plt.title("miRNA–Clinical Correlations   [" + " • ".join(title_bits) + "]", fontsize=14)
    plt.xlabel("Clinical Variables", fontsize=12)
    plt.ylabel("miRNAs", fontsize=12)
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)

    # Save
    plt.tight_layout()
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    print(f"Heatmap saved to: {args.output}")


if __name__ == "__main__":
    main()
