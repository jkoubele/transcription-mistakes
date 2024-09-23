import json
import pandas as pd
from pathlib import Path
from tqdm import tqdm
from collections import defaultdict
from typing import NamedTuple



class Observation(NamedTuple):
    base: str
    base_next: str
    base_prev: str
    phase: str
    gsm: str

cell_cycle_info_folder = Path('/cellfile/datapublic/jkoubele/sc_data/cell_cycle_info')
detected_errors_folder = Path('/cellfile/datapublic/jkoubele/sc_data/detected_errors_3')

gse = 'GSE190848'
# gsm = 'GSM4425803'

rows_correct = defaultdict(int)
rows_mismatch = defaultdict(int)
gse_folder = detected_errors_folder / gse
for gsm_folder in gse_folder.iterdir():
    print(gsm_folder)
    gsm = gsm_folder.stem
    
    metadata = pd.read_csv(cell_cycle_info_folder / gse / gsm / 'metadata.tsv', sep='\t')
    metadata = metadata[metadata['percent.mt']<5.0]
    metadata = metadata.set_index('rowname', drop=False)
    gsm_json_path = detected_errors_folder / gse / gsm
    for subfolder in gsm_json_path.iterdir():
        print(subfolder)
        for file_path in tqdm(subfolder.iterdir()):
            with open(file_path) as file:
                detected_errors_json = json.load(file)
            cell_barcode = file_path.stem
            if cell_barcode not in metadata.index:
                continue
            for features, base_count in detected_errors_json:            
                observation = Observation(base=features['base'],
                                          base_next=features['base_next'],
                                          base_prev=features['base_prev'],
                                          phase=metadata.loc[cell_barcode]['Phase'],
                                          gsm=gsm)
                
                num_correct = base_count[features['base']]
                num_mismatch = sum([value for key, value in base_count.items() if key != features['base']])
                
                if num_correct > 0:
                    rows_correct[observation] += num_correct
                if num_mismatch > 0:
                    rows_mismatch[observation] += num_mismatch                
        
    
rows = []
for key, value in rows_correct.items():
    row = key._asdict()
    row['y'] = 0
    row['count'] = value
    rows.append(row)
for key, value in rows_mismatch.items():
    row = key._asdict()
    row['y'] = 1
    row['count'] = value
    rows.append(row)
    
df = pd.DataFrame.from_records(rows) 

output_folder = Path('/cellfile/datapublic/jkoubele/sc_data/design_matrices') / gse
output_folder.mkdir(exist_ok=True, parents=True)
df.to_csv(output_folder / 'design_matrix.tsv', sep='\t', index=False)


