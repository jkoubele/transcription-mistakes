import os

os.environ[
    'OPENBLAS_NUM_THREADS'] = '1'  # solves weird error when importing numpy (and consequently e.g. pandas, biopython etc.) on cluster

import argparse
from collections import defaultdict
from dataclasses import dataclass, field
from itertools import chain
from typing import Optional
import pysam
from Bio import SeqIO, Seq
from dataclasses_json import DataClassJsonMixin
from tqdm import tqdm
from pathlib import Path
from time import time

MIN_LENGTH_OF_VALID_CB_OR_UB_TAG = 5

DNA_ALPHABET = ('A', 'C', 'G', 'T')

@dataclass
class ReadWithPosition:
    read: pysam.AlignedSegment
    position: int
    
    

if __name__ == "__main__":    
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_folder',
                        default='/cellfile/datapublic/jkoubele/sc_data/split_by_cell/GSE147319/GSM4425801/AA',
                        help='Folder containing .bam files to be processed.')
    parser.add_argument('--output_folder',
                        default='/cellfile/datapublic/jkoubele/sc_data/detected_errors/GSE147319/GSM4425801',
                        help='Output JSON file will be written here.')
    parser.add_argument('--reference_genome_fasta_file',
                        default='/cellfile/datapublic/jkoubele/reference_genomes/GRCm39/Mus_musculus.GRCm39.dna.primary_assembly.fa',                       
                        help='FASTA file of the reference genome that was used for the alignment of the input .bam files.')
    args = parser.parse_args()
    Path(args.output_folder).mkdir(exist_ok=True, parents=True)
    input_folder = Path(args.input_folder)
    #%%
    reference_genome_all_sequences = list(SeqIO.parse(open(args.reference_genome_fasta_file), 'fasta'))
    reference_genome: dict[str, Seq] = {record.name: record.seq for record in reference_genome_all_sequences
                                        if record.name in ('X', 'Y') or record.name.isdigit()}
    
    
    #%%    
    def find_consensus_base(reads: list[ReadWithPosition]) -> Optional[str]:
        unique_bases = set([x.read.query_sequence[x.position] for x in reads])
        if len(unique_bases) != 1:
            return None
        base = unique_bases.pop()
        return base if base in DNA_ALPHABET else None
    
    
    min_reads_per_evaluated_umi = 3
    min_umis_to_check_for_mutation = 10
    for file in tqdm([file for file in input_folder.iterdir() if file.suffix == '.bam']):
        num_evaluated = 0        
        num_evaluated_valid_base = 0
        num_mismatches = 0
        time_start = time()
        bamfile = pysam.AlignmentFile(file, "rb")
        
        # seen_umis = set()
        # for read in bamfile:
        #     if read.mapping_quality < 255 or not read.has_tag('UB'):
        #         continue   
        #     umi = read.get_tag('UB')  
        #     if umi in seen_umis:
        #         continue
        #     seen_umis.add(umi)
        #     num_evaluated += len([x for x in read.get_aligned_pairs() if x[0] is not None and x[1] is not None])
           
        # break
        
        for chromosome_name in tqdm(reference_genome.keys()):
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
                    for umi, reads in reads_by_umi.items():
                        if len(reads) >= min_reads_per_evaluated_umi:
                            # evaluate for polymerase error
                            num_evaluated += 1
                            consensus_base = find_consensus_base(reads)
                            if consensus_base is not None:
                                num_evaluated_valid_base += 1
                                if consensus_base != reference_symbol:
                                    num_mismatches += 1
        del pileup_column
        # del pileup_read
        time_end = time()
        print()
        print(f"{num_evaluated / 1e3}k")
        print(f"Time: {time_end - time_start}")
        break
    