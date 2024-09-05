import requests
from time import sleep
import xml.etree.ElementTree as ET
import json
from tqdm import tqdm
import pandas as pd
from pathlib import Path


def get_response(url):
    sleep(0.1)
    return requests.get(url)


if __name__ == "__main__":
    einfo_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi?retmode=json"
    einfo_json_sra = get_response(einfo_url + "&db=sra").json()
    einfo_json_gds = get_response(einfo_url + "&db=gds").json()

    sc_datasets = pd.read_csv(Path('metadata/sc_datasets.tsv'), sep='\t')
    gse_list = sc_datasets['GSE'].to_list()

    res = {}

    for gse in tqdm(gse_list, desc='Querying eutils'):

        search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?api_key=bcad3816f864fe73bd100b1c8fb72b349008&retmode=json&"
        search_query = f"db=gds&term={gse}[GEO Accession]+AND+gse[Entry Type]"
        search_response_json = get_response(search_url + search_query).json()
        returned_id_list = search_response_json['esearchresult']['idlist']

        assert len(returned_id_list) == 1
        gse_uid = returned_id_list[0]

        samples_info = []

        summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?api_key=bcad3816f864fe73bd100b1c8fb72b349008&retmode=json&"
        summary_query = f"db=gds&id={gse_uid}"
        summary_response_json = get_response(summary_url + summary_query).json()
        gse_summary = summary_response_json['result'][gse_uid]
        for sample in gse_summary['samples']:
            gsm = sample['accession']
            search_query_sample = f"db=gds&term={gsm}[GEO Accession]+AND+gsm[Entry Type]"
            search_response_sample_json = get_response(search_url + search_query_sample).json()
            sample_id_list = search_response_sample_json['esearchresult']['idlist']
            assert len(sample_id_list) == 1

            gsm_uid = sample_id_list[0]
            summary_sample_query = f"db=gds&id={gsm_uid}"
            summary_sample_response_json = get_response(summary_url + summary_sample_query).json()
            gsm_summary = summary_sample_response_json['result'][gsm_uid]
            taxon = gsm_summary['taxon']
            sra_relations = [relation for relation in gsm_summary['extrelations'] if relation['relationtype'] == 'SRA']
            assert len(sra_relations) == 1

            srx = sra_relations[0]['targetobject']
            assert srx.startswith('SRX')

            search_query_srx = f"db=sra&term={srx}[Accession]"
            search_response_srx_json = get_response(search_url + search_query_srx).json()
            srx_id_list = search_response_srx_json['esearchresult']['idlist']
            assert len(srx_id_list) == 1
            srx_uid = srx_id_list[0]

            summary_srx_query = f"db=sra&id={srx_uid}"
            summary_srx_response_json = get_response(summary_url + summary_srx_query).json()
            srx_summary = summary_srx_response_json['result'][srx_uid]

            runs_xml = ET.fromstring(f"<Runs>{srx_summary['runs'].strip()}</Runs>")
            for run in runs_xml:
                srr = run.attrib['acc']
                assert srr.startswith('SRR')
                is_public = run.attrib['is_public']

            samples_info.append({'GSM': gsm,
                                 'SRX': srx,
                                 'taxon': taxon,
                                 'title': gsm_summary['title'],
                                 'summary': gsm_summary['summary'],
                                 'runs_srr': [run.attrib['acc'] for run in runs_xml]})
        res[gse] = {'title': gse_summary['title'],
                    'summary': gse_summary['summary'],
                    'samples': samples_info}

    with open(Path('metadata/datasets_overview.json'), 'w') as file:
        json.dump(res, file)
