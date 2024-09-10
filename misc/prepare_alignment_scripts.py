import json
from pathlib import Path

from paths import repository_path

if __name__ == "__main__":
    with open(repository_path / Path('misc/alignment_script_template.sh')) as file:
        alignment_script_template = file.read()

    with open(repository_path / Path('metadata/datasets_overview.json')) as file:
        datasets_overview = json.load(file)

    barcodes_file_3m = '3M-february-2018.txt'
    barcodes_file_737k = '737k-august-2016.txt'
    for gse, barcodes_file, umi_length in [('GSE190848', barcodes_file_3m, str(12)),
                                           ('GSE230816', barcodes_file_3m, str(12)),
                                           ('GSE147319', barcodes_file_737k, str(10)),
                                           ('GSE158859', barcodes_file_3m, str(12))]:
        dataset = datasets_overview[gse]
        samples = {sample['GSM']: sample for sample in dataset['samples']}

        for gsm, sample in samples.items():
            reads_1 = ','.join([f"/sra_folder/{srr_id}/{srr_id}_1.fastq.gz" for srr_id in sample['runs_srr']])
            reads_2 = ','.join([f"/sra_folder/{srr_id}/{srr_id}_2.fastq.gz" for srr_id in sample['runs_srr']])

            alignment_script = alignment_script_template
            for placeholder, value in [('GSE', gse),
                                       ('GSM', gsm),
                                       ('BARCODES_FILE', barcodes_file),
                                       ('READS_1', reads_1),
                                       ('READS_2', reads_2),
                                       ('UMI_LENGTH', umi_length)]:
                alignment_script = alignment_script.replace(placeholder, value)
            output_dir = repository_path / Path(f'alignment_scripts/{gse}/{gsm}')
            output_dir.mkdir(parents=True, exist_ok=True)
            with open(output_dir / 'alignment_script.sh', 'w') as file:
                file.write(alignment_script)
