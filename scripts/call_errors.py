import os

os.environ[
    'OPENBLAS_NUM_THREADS'] = '1'  # solves weird error when importing numpy (and consequently e.g. pandas, biopython etc.) on cluster

import argparse
from collections import defaultdict
from dataclasses import dataclass
from typing import Optional, NamedTuple
import pysam
from Bio import SeqIO, Seq
from tqdm import tqdm
from pathlib import Path
from time import time
import json

MIN_LENGTH_OF_VALID_CB_OR_UB_TAG = 5

DNA_ALPHABET = ('A', 'C', 'G', 'T')


@dataclass
class ReadWithPosition:
    read: pysam.AlignedSegment
    position: int


class Features(NamedTuple):
    base: str
    base_next: str
    base_prev: str


def convert_counter_to_serializable(counter: dict) -> list:
    return [[key._asdict(), value]
            for key, value in counter.items()]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_folder',
                        default='/cellfile/datapublic/jkoubele/sc_data/split_by_cell/GSE147319/GSM4425801/AA',
                        help='Folder containing .bam files to be processed.')
    parser.add_argument('--output_folder',
                        default='/cellfile/datapublic/jkoubele/sc_data/detected_errors/GSE147319/GSM4425801',
                        help='Output JSON file will be written here.')
    parser.add_argument('--output_folder_locations',
                        default='/cellfile/datapublic/jkoubele/sc_data/mismatch_locations/GSE147319/GSM4425801',
                        help='Output JSON file will be written here.')
    parser.add_argument('--reference_genome_fasta_file',
                        default='/cellfile/datapublic/jkoubele/reference_genomes/GRCm39/Mus_musculus.GRCm39.dna.primary_assembly.fa',
                        help='FASTA file of the reference genome that was used for the alignment of the input .bam files.')
    args = parser.parse_args()
    
    input_folder = Path(args.input_folder)
    
    output_folder = Path(args.output_folder)
    output_folder.mkdir(exist_ok=True, parents=True)
    
    output_folder_locations = Path(args.output_folder_locations)
    output_folder_locations.mkdir(exist_ok=True, parents=True)
    # %%
    reference_genome_all_sequences = list(SeqIO.parse(open(args.reference_genome_fasta_file), 'fasta'))
    reference_genome: dict[str, Seq] = {record.name: record.seq for record in reference_genome_all_sequences
                                        if record.name in ('X', 'Y') or record.name.isdigit()}


    # %%
    def find_consensus_base(reads: list[ReadWithPosition]) -> Optional[str]:
        unique_bases = set([x.read.query_sequence[x.position] for x in reads])
        if len(unique_bases) != 1:
            return None
        base = unique_bases.pop()
        return base if base in DNA_ALPHABET else None


    min_reads_per_evaluated_umi = 5
    min_umis_to_check_for_mutation = 10
    for file in tqdm([file for file in input_folder.iterdir() if file.suffix == '.bam']):        
        time_start = time()

        counter = defaultdict(lambda: {base: 0 for base in DNA_ALPHABET})
        mismatch_locations = []

        bamfile = pysam.AlignmentFile(file, "rb")

        for chromosome_name in tqdm(reference_genome.keys(), disable=True):
            for pileup_column in bamfile.pileup(chromosome_name):
                if pileup_column.get_num_aligned() < min_reads_per_evaluated_umi + min_umis_to_check_for_mutation:
                    continue
                forward_reads_by_umi: dict[str, list[ReadWithPosition]] = defaultdict(list)
                reverse_reads_by_umi: dict[str, list[ReadWithPosition]] = defaultdict(list)

                reference_symbol = reference_genome[chromosome_name][pileup_column.reference_pos]

                for pileup_read in pileup_column.pileups:
                    if pileup_read.query_position is not None:
                        read = pileup_read.alignment

                        if read.mapping_quality < 255 or not read.has_tag('UB'):
                            continue
                        umi = read.get_tag('UB')
                        if len(umi) < MIN_LENGTH_OF_VALID_CB_OR_UB_TAG:
                            continue
                        if read.is_forward:
                            forward_reads_by_umi[umi].append(ReadWithPosition(read=read,
                                                                              position=pileup_read.query_position))
                        else:
                            reverse_reads_by_umi[umi].append(ReadWithPosition(read=read,
                                                                              position=pileup_read.query_position))
                for is_forward_strand, reads_by_umi in [(True, forward_reads_by_umi),
                                                        (False, reverse_reads_by_umi)]:
                    if len(reads_by_umi) < min_umis_to_check_for_mutation:
                        continue
                    high_coverage_umi = {key: value for key, value in reads_by_umi.items() if
                                         len(value) >= min_reads_per_evaluated_umi}
                    low_coverage_umi = {key: value for key, value in reads_by_umi.items() if
                                        len(value) < min_reads_per_evaluated_umi}

                    consensus_bases_high_coverage = [find_consensus_base(reads) for reads in high_coverage_umi.values()]
                    consensus_bases_high_coverage = [base for base in consensus_bases_high_coverage if base is not None]

                    mismatches_high_coverage = [base for base in consensus_bases_high_coverage if
                                                base != reference_symbol]

                    if not mismatches_high_coverage:
                        base_next = reference_genome[chromosome_name][
                            pileup_column.reference_pos + 1] if is_forward_strand else \
                            reference_genome[chromosome_name][pileup_column.reference_pos - 1]
                        base_prev = reference_genome[chromosome_name][
                            pileup_column.reference_pos - 1] if is_forward_strand else \
                            reference_genome[chromosome_name][pileup_column.reference_pos + 1]
                        features = Features(base=reference_symbol,
                                            base_next=base_next,
                                            base_prev=base_prev)
                        counter[features][reference_symbol] += len(consensus_bases_high_coverage)

                    elif len(mismatches_high_coverage) == 1:
                        mismatch_base = mismatches_high_coverage[0]
                        # check for evidence of somatic mutation
                        consensus_bases_low_coverage = [find_consensus_base(reads) for reads in
                                                        low_coverage_umi.values()]
                        if not any([base == mismatch_base for base in consensus_bases_low_coverage]):
                            base_next = reference_genome[chromosome_name][
                                pileup_column.reference_pos + 1] if is_forward_strand else \
                                reference_genome[chromosome_name][pileup_column.reference_pos - 1]
                            base_prev = reference_genome[chromosome_name][
                                pileup_column.reference_pos - 1] if is_forward_strand else \
                                reference_genome[chromosome_name][pileup_column.reference_pos + 1]
                            features = Features(base=reference_symbol,
                                                base_next=base_next,
                                                base_prev=base_prev)
                            counter[features][mismatch_base] += 1
                            mismatch_locations.append({'chr': chromosome_name,
                                                       'position': pileup_column.reference_pos,
                                                       'strand': '+' if is_forward_strand else '-',
                                                       'reference_base': reference_symbol,
                                                       'mismatch_base': mismatch_base,
                                                       'num_umi': len(reads_by_umi)})

                        pass
                    else:
                        pass  # too many mismatches

                    consensus_bases_high_coverage = [find_consensus_base(reads) for reads in high_coverage_umi.values()]
                    mismatches_high_coverage = [base for base in consensus_bases_high_coverage if
                                                base != reference_symbol]

                    # for umi, reads in reads_by_umi.items():
                    #     if len(reads) >= min_reads_per_evaluated_umi:
                    #         # evaluate for polymerase error
                    #         num_evaluated += 1
                    #         consensus_base = find_consensus_base(reads)
                    #         if consensus_base is not None:
                    #             num_evaluated_valid_base += 1
                    #             if consensus_base != reference_symbol:
                    #                 num_mismatches += 1
                    # base_next = reference_genome[chromosome_name][
                    #     pileup_column.reference_pos + 1] if is_forward_strand else \
                    #     reference_genome[chromosome_name][pileup_column.reference_pos - 1]
                    # base_prev = reference_genome[chromosome_name][
                    #     pileup_column.reference_pos - 1] if is_forward_strand else \
                    #     reference_genome[chromosome_name][pileup_column.reference_pos + 1]
                    # features = Features(base=reference_symbol,
                    #                     base_next=base_next,
                    #                     base_prev=base_prev)
                    # counter[features][consensus_base] += 1

        del pileup_column

        # del pileup_read

        with open(output_folder / f"{file.stem}.json", 'w') as file_out:
            json.dump(convert_counter_to_serializable(counter), file_out)

        if mismatch_locations:
            with open(output_folder_locations / f"{file.stem}.json", 'w') as file_out:
                json.dump(mismatch_locations, file_out)

        time_end = time()
        print()        
        print(f"Time: {time_end - time_start}")
