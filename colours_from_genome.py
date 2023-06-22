from genome_functions import read_pombe_genome, make_synonym_dict
from custom_biopython import SeqIO
import pandas
import re
import glob
import os
import sys

synonym_dict = make_synonym_dict('valid_ids_data/gene_IDs_names.tsv', 'valid_ids_data/obsoleted_ids.tsv', 'valid_ids_data/missing_synonyms.tsv')

def count_colours(revision_folder):
    counts_file = f'{revision_folder}/counts.txt'
    if os.path.isfile(counts_file):
        return
    colours = list()
    for contig_file in glob.glob(f'{revision_folder}/*.contig'):
        print(f'    {contig_file}')
        genome = read_pombe_genome(contig_file, 'embl', synonym_dict, 'valid_ids_data/all_systematic_ids_ever.txt','valid_ids_data/known_exceptions.tsv')
        cds_features = [feature for feature in genome.features if feature.type == 'CDS']

        locus_ids = set(pandas.read_csv('valid_ids_data/gene_IDs_names.tsv', delimiter='\t', na_filter=False, dtype=str)['systematic_id'])

        for feature in cds_features:
            if 'systematic_id' in feature.qualifiers:
                systematic_id = feature.qualifiers['systematic_id'][0]
            else:
                systematic_id = None
                print('Feature without systematic id')
                print(feature)
            # If multi-transcript, skip if not ending in .1
            if systematic_id is not None and systematic_id not in locus_ids and re.sub(r'\.\d$', '', systematic_id) in locus_ids and not systematic_id.endswith('.1'):
                continue
            if 'colour' not in feature.qualifiers:
                print(f'colour missing in {systematic_id}')
                continue
            # Sometimes there are repeated colours
            colour = list(set(feature.qualifiers['colour']))
            if len(colour) != 1:
                print(f'two colours in {systematic_id}')
                # The minimal value is the highest level of characterisation
                colour = [min(colour)]
            colours.append(int(colour[0]))


    with open(counts_file, 'w') as out:
        out.write(' '.join(str(colours.count(i)) for i in range(16)))
        out.write('\n')


if __name__ == '__main__':
    for folder in sys.argv[1:]:
        print(f'Counting for {folder}')
        count_colours(folder)
