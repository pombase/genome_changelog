from Bio import Entrez
from custom_biopython import SeqIO
import os
import requests
import pandas

data = pandas.read_csv('results/genome_sequence_changes.tsv', sep='\t')
data.fillna('', inplace=True)
data = data[data['reference'] != ''].copy()

if not os.path.exists('test_folder'):
    os.mkdir('test_folder')
if not os.path.exists('test_folder/accession_numbers/'):
    os.mkdir('test_folder/accession_numbers/')

for accession, link in zip(data['reference'], data['link']):

    print(f'\033[0;32mComparing {accession} ========================================\033[0m')

    # Download the genbank file from NCBI
    handle = Entrez.efetch(db='nucleotide', id=accession, rettype='gb', retmode='text')
    filename = f'test_folder/accession_numbers/{accession}.gb'
    with open(filename, 'w') as outfile:
        outfile.write(handle.read())

    # Download the corresponding one from pombase

    with open(f'test_folder/accession_numbers/{accession}.embl', 'w') as out:
        out.write(requests.get(link).text)

    # Verify that their sequences are identical

    pombase_file = SeqIO.read(f'test_folder/accession_numbers/{accession}.embl', 'embl')
    genbank_file = SeqIO.read(f'test_folder/accession_numbers/{accession}.gb', 'gb')

    assert(pombase_file.seq == genbank_file.seq)



