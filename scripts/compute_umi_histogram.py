import json
import argparse
from collections import defaultdict, Counter
import pysam
from pathlib import Path
import pandas as pd
from tqdm import tqdm
import numpy as np

MIN_VALID_UMI_LENGTH = 8

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_folder')
    parser.add_argument('--output_folder')
    args = parser.parse_args()

    input_folder = Path(args.input_folder)
    output_folder = Path(args.output_folder)

    output_folder.mkdir(parents=True, exist_ok=True)

    barcodes_df = pd.read_csv(input_folder / Path('Solo.out/GeneFull/filtered/barcodes.tsv'),
                              sep='\t',
                              names=['barcode'])

    umis_by_cell = {cell_barcode: defaultdict(int) for cell_barcode in barcodes_df['barcode']}

    bamfile = pysam.AlignmentFile(input_folder / 'Aligned.sortedByCoord.out.bam', "rb")
    for read in tqdm(bamfile):
        cell_barcode = read.get_tag('CB')
        umi = read.get_tag('UB')
        if len(umi) < MIN_VALID_UMI_LENGTH:
            continue  # filtering out reads with empty UMI tag, which can be '-' or similar
        if cell_barcode in umis_by_cell:
            umis_by_cell[cell_barcode][umi] += 1

    # %%
    output_overall = {}
    output_by_cell = {}

    histograms = {barcode: dict(sorted(Counter(umi_counts.values()).items()))
                  for barcode, umi_counts in umis_by_cell.items()}
    cumulative_histograms = {barcode: {key: value for key, value in
                                       zip(histogram.keys(),
                                           np.flip(np.cumsum(list(reversed((histogram.values()))))).tolist())} for
                             barcode, histogram in histograms.items()}

    overall_histogram = defaultdict(int)
    overall_cumulative_histogram = defaultdict(int)

    for histogram in histograms.values():
        for key, value in histogram.items():
            overall_histogram[key] += value

    for histogram in cumulative_histograms.values():
        for key, value in histogram.items():
            overall_cumulative_histogram[key] += value

    output_overall['overall_histogram'] = dict(sorted(overall_histogram.items()))
    output_overall['overall_cumulative_histogram'] = dict(sorted(overall_cumulative_histogram.items()))

    output_by_cell['histograms'] = histograms
    output_by_cell['cumulative_histograms'] = cumulative_histograms

    with open(output_folder / 'umi_histograms_overall.json', 'w') as file:
        json.dump(output_overall, file)

    with open(output_folder / 'umi_histograms_by_cell.json', 'w') as file:
        json.dump(output_by_cell, file)
