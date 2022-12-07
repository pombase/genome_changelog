"""
Make a dictionary tsv which links obsoleted identifiers to the current identifiers
"""

import glob
from genome_functions import read_pombe_genome, make_synonym_dict

print('obsolete_id','systematic_id',sep='\t')

for contig_file in 'chromosome1 chromosome2 chromosome3 mating_type_region pMIT'.split():

    # Get all revisions downloaded and use the latest one to read the info
    last_file = max(glob.glob(f'data/{contig_file}/*.contig'),key=lambda x: int(x.split('/')[-1].split('.')[0]))
    synonym_dictionary = make_synonym_dict('valid_ids_data/gene_IDs_names.tsv')
    contig = read_pombe_genome(last_file,'embl',synonym_dictionary, 'valid_ids_data/all_systematic_ids_ever.txt','valid_ids_data/known_exceptions.tsv')

    for feature in contig.features:
        if 'obsolete_name' in feature.qualifiers:
            for obsolete_id in feature.qualifiers['obsolete_name']:
                if obsolete_id.startswith('SP'):
                    print(obsolete_id, feature.qualifiers['systematic_id'][0],sep='\t')
