#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
DESCRIPTION:
Filters out miRNAs with low expression based on median counts across groups.
Groups are defined dynamically using a status CSV file (sample,status).
Only miRNAs with a median ≥ 10 in at least one group are retained.

Also identifies miRNAs with low variability (fewer than 20 unique values).

INPUT:
- A tab-delimited count matrix file (first col: GeneID, rest: sample columns)
- A CSV file mapping sample names to group status
OUTPUT:
- A file with retained miRNAs and their group-wise medians
- A file with retained miRNAs in original count matrix format

USAGE:
python filter_out_low_counts.py countData.txt lipid_status.csv kept_medians.txt kept_data.txt
"""

import sys
import statistics
import csv

def load_status_map(status_file):
    """Load sample-to-status mapping from CSV."""
    status_map = {}
    with open(status_file, encoding="utf-8-sig") as f:  # utf-8-sig removes BOM if present
        reader = csv.DictReader(f)
        # Normalize fieldnames (strip spaces, lowercase)
        reader.fieldnames = [fn.strip().lower() for fn in reader.fieldnames]
        for row in reader:
            sample = row["sample"].strip()
            status = row["status"].strip()
            status_map[sample] = status
    return status_map


def main():
    if len(sys.argv) != 5:
        print("USAGE: python filter_out_low_counts.py INPUT_MATRIX STATUS_FILE OUTPUT_MEDIANS OUTPUT_MATRIX")
        sys.exit(1)

    rawCountFileName = sys.argv[1]
    statusFile = sys.argv[2]
    mediansOutputPath = sys.argv[3]
    matrixOutputPath = sys.argv[4]
    low_variability_threshold = 1  # fewer than 20 unique values = low variability

    # Load status mapping
    status_map = load_status_map(statusFile)

    try:
        with open(rawCountFileName) as rawCounts, \
             open(mediansOutputPath, 'w') as keptWithMediansFile, \
             open(matrixOutputPath, 'w') as keptDataFile:

            # --- Read header ---
            header = rawCounts.readline().strip().split("\t")
            sample_names = header[1:]  # skip GeneID

            # --- Build group indices dynamically ---
            groups = {}
            for i, sample in enumerate(sample_names, start=1):
                if sample in status_map:
                    group = status_map[sample]
                    groups.setdefault(group, []).append(i)
                else:
                    print(f"Warning: Sample '{sample}' not found in status file")

            # --- Write header lines ---
            group_labels = sorted(groups.keys())
            medianHeader = "GeneID\t" + "\t".join(f"MedianGroup{g}" for g in group_labels) + "\n"
            keptWithMediansFile.write(medianHeader)
            keptDataFile.write("\t".join(header) + '\n')

            # --- Process each miRNA row ---
            for line in rawCounts:
                fields = line.strip().split("\t")
                geneId = fields[0]
                values = [float(x) for x in fields[1:]]

                # Skip low variability
                if len(set(values)) < low_variability_threshold:
                    continue

                # Compute medians per group
                group_medians = []
                for g in group_labels:
                    indices = groups[g]
                    group_counts = [values[i-1] for i in indices]  # adjust index
                    group_median = statistics.median(group_counts)
                    group_medians.append(group_median)

                # Keep if median ≥ 10 in any group
                if any(m >= 10 for m in group_medians):
                    medianLine = f'{geneId}\t' + "\t".join(map(str, group_medians)) + '\n'
                    keptWithMediansFile.write(medianLine)
                    keptDataFile.write(line)

    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(2)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        pass
