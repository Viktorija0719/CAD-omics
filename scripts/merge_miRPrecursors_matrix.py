#!/usr/bin/env python3

"""
DESCRIPTION:
This script processes miRNA expression files from mirDeep2, merges counts from different precursors 
of the same mature miRNA across multiple samples, then converts the merged data into a matrix format 
(Sample_Id as columns, miRNAs as rows). 

It also formats numeric values: integers are saved without ".0" while ".5" values are preserved.

USAGE:
    python merge_miRPrecursors_matrix.py <input_dir> <sample_id_file.csv> <output_file.txt>

ARGUMENTS:
    <input_dir>          Directory containing files like 'miRNAs_expressed_all_samples_*.csv'
    <sample_id_file.csv> CSV mapping CeGaT IDs to clinical sample IDs
    <output_file.txt>    Output file path for the final matrix
"""

import os
import sys
import glob
import statistics
import pandas as pd


def main():
    # --- Check arguments ---
    if len(sys.argv) < 4:
        sys.exit(f"USAGE: {sys.argv[0]} <miRNA_input_dir> <sample_id_file.csv> <output_file_path>")

    miRNA_dir, sample_id_file, output_file_path = sys.argv[1:4]

    # --- Load clinical ID map ---
    clinIDs = load_clinical_ids(sample_id_file)

    # --- Find input files ---
    allFiles = glob.glob(os.path.join(miRNA_dir, 'miRNAs_expressed_all_samples_*.csv'))
    if not allFiles:
        sys.exit(f"No input files found in '{miRNA_dir}'")

    # --- Step 1: Merge precursor counts ---
    counts_map = {}
    for file in sorted(allFiles):
        sample = extract_sample_name(file)
        with open(file) as f:
            next(f)  # skip header
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) == 6:
                    miRNA, count, precursor = parts[0], parts[1], parts[2]
                    counts_map.setdefault(f"{miRNA}|{precursor}", []).append((sample, float(count)))
                else:
                    print(f"Warning: skipped malformed line in {file}: {line.strip()}")

    # Condense by averaging counts for same miRNA+sample
    condensed = {}
    for key, values in counts_map.items():
        miRNA = key.split("|")[0]
        for sample, count in values:
            condensed.setdefault(f"{miRNA}|{sample}", []).append(count)

    merged_data = []
    for key in sorted(condensed):
        miRNA, sample = key.split("|")
        if sample in clinIDs:
            mean_count = statistics.mean(condensed[key])
            merged_data.append([miRNA, clinIDs[sample], mean_count])

    # --- Step 2: Convert to DataFrame ---
    df = pd.DataFrame(merged_data, columns=["miRNA", "Sample_Id", "Count"])

    # --- Step 3: Pivot to matrix format ---
    df = df.sort_values(by="Sample_Id")
    mat = df.pivot(index="miRNA", columns="Sample_Id", values="Count").reset_index()

    # --- Step 4: Format numeric values ---
    mat = mat.apply(pd.to_numeric, errors='ignore')

    def format_value(val):
        if isinstance(val, float):
            if val.is_integer():
                return str(int(val))       # remove .0
            elif val * 10 % 5 == 0:        # keep .5 values
                return str(val)
            else:
                return str(val)
        return val

    mat = mat.applymap(format_value)

    # --- Step 5: Save final table ---
    mat.to_csv(output_file_path, sep="\t", index=False)
    print(f"✅ Output saved to {output_file_path}")


def extract_sample_name(filepath):
    """Extract CeGaT sample name from filename."""
    base = os.path.basename(filepath)
    return base.split(".")[0].split("_")[-1].replace("miRNA", "")


def load_clinical_ids(path):
    """Load CeGaT → Clinical ID mapping from CSV."""
    id_map = {}
    with open(path) as f:
        next(f)  # skip header
        for line in f:
            cegat_id, clin_id = line.strip().split(",")[:2]
            id_map[cegat_id] = clin_id
    return id_map


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        pass
