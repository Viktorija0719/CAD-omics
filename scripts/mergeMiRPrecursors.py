#!/usr/bin/env python3


# DESCRIPTION:
# This script processes miRNA expression files produced by mirDeep2, combining counts 
# from different precursors of the same mature miRNA across multiple samples.
# It outputs a table where each row is a miRNA-clinical sample pair with the averaged read count.
#
# USAGE:
#   python 05_mergeMiRPrecursors.py <input_dir> <sample_id_file.csv> <output_file>
#
# ARGUMENTS:
#   <input_dir>          Directory containing files named like 'miRNAs_expressed_all_samples_*.csv'
#   <sample_id_file.csv> CSV file mapping internal sample IDs (CeGaT) to clinical sample IDs
#   <output_file>        Path to the output file (e.g., ./countTable_11.txt)
# NOTE:
#   This script was adapted and modified from original code developed by 
#   Dr. rer. nat. Baiba Vilne.


import os, sys, glob, statistics

def main():
    if len(sys.argv) < 4:
        print(f"USAGE: {sys.argv[0]} <miRNA_input_dir> <sample_id_file.csv> <output_file_path>")
        sys.exit(1)

    miRNA_dir, sample_id_file, output_file_path = sys.argv[1:4]
    clinIDs = load_clinical_ids(sample_id_file)
    allFiles = glob.glob(os.path.join(miRNA_dir, 'miRNAs_expressed_all_samples_*.csv'))

    if not allFiles:
        print(f"️  No input files found in '{miRNA_dir}'")
        sys.exit(1)

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
                    print(f"Warning: skipped line in {file}: {line.strip()}")

    condensed = {}
    for key, values in counts_map.items():
        miRNA = key.split("|")[0]
        for sample, count in values:
            condensed.setdefault(f"{miRNA}|{sample}", []).append(count)

    with open(output_file_path, "w") as out:
        print("miRNA\tSample_Id\tCount", file=out)
        for key in sorted(condensed):
            miRNA, sample = key.split("|")
            if sample in clinIDs:
                mean_count = statistics.mean(condensed[key])
                print(f"{miRNA}\t{clinIDs[sample]}\t{mean_count}", file=out)

def extract_sample_name(filepath):
    # Pulls last underscore-separated part before `.csv`
    base = os.path.basename(filepath)
    return base.split(".")[0].split("_")[-1].replace("miRNA", "")

def load_clinical_ids(path):
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