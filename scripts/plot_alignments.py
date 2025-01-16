#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt

def main(input_path, output_dir):
	
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)

	alnSummary = pd.read_csv(input_path, sep="\t", skiprows=[1])

	# plot alignment positions colored by idy
	plt.figure(figsize=(16, 8))
	chromosomes = alnSummary['#chr'].unique()

	# make y chromosome positions sorted by chr number
	chromosomes = sorted(alnSummary['#chr'].unique(), key=lambda x: (int(x[3:]) if x[3:].isdigit() else float('inf')))
	chrom_map = {chrom: i+1 for i, chrom in enumerate(chromosomes)}

	# plot lines for each position
	for _, row in alnSummary.iterrows():
	    plt.plot([row['start_pos'], row['end_pos']],
	             [chrom_map[row['#chr']], chrom_map[row['#chr']]], 
	             color=plt.cm.viridis(row['identity']), linewidth=10, alpha=1)

	# labels
	plt.gca().invert_yaxis()
	plt.yticks(ticks=list(chrom_map.values()), labels=list(chrom_map.keys()))
	plt.xlabel("Position")
	plt.ylabel("Chromosome")
	plt.title("BED File Positions by Chromosome and Identity")
	plt.colorbar(plt.cm.ScalarMappable(cmap='viridis'), label="Identity")

	plt.tight_layout()
	plt.savefig(output_dir+"/alignment_summary.png", dpi=300)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-a","--alignment_summary",
        required=True,
        type=str,
        help="Input file of alignment_summary.tsv"
    )

    parser.add_argument(
        "-o","--output_dir",
        required=True,
        type=str,
        help="Directory path where output will be written"
    )

    args = parser.parse_args()

    main(input_path=args.alignment_summary, output_dir=args.output_dir)

    