from genome_functions import read_pombe_genome, make_synonym_dict
from custom_biopython import SeqIO
import pandas
import re
# We override this method to allow no space between number and BP

synonym_dict = make_synonym_dict('valid_ids_data/gene_IDs_names.tsv', 'valid_ids_data/obsoleted_ids.tsv', 'valid_ids_data/missing_synonyms.tsv')

def parse_genome(genome_file, date):
    # Special reader for old genomes handles several edge-cases
    if date < '2023-03-01':
        return read_pombe_genome(genome_file, 'embl', synonym_dict, 'valid_ids_data/all_systematic_ids_ever.txt','valid_ids_data/known_exceptions.tsv')
    else:
        with open(genome_file, errors='replace') as ins:
            return SeqIO.read(ins,'embl')


colours = list()

for contig_file in ['data/chromosome1/8752.contig']:
    genome = parse_genome(contig_file, '2009-07-30')
    cds_features = [feature for feature in genome.features if feature.type == 'CDS']

    locus_ids = set(pandas.read_csv('valid_ids_data/gene_IDs_names.tsv', delimiter='\t', na_filter=False, dtype=str)['systematic_id'])

    for feature in cds_features:
        # If multi-transcript, skip if not ending in .1
        systematic_id = feature.qualifiers['systematic_id'][0]
        if systematic_id not in locus_ids and re.sub(r'\.\d$', '', systematic_id) in locus_ids and not systematic_id.endswith('.1'):
            continue

        colour = feature.qualifiers['colour']
        if len(colour) != 1:
            raise ValueError(f'two colours in {systematic_id}')
        colours.append(int(colour[0]))


with open('gene_characterisation/counts.txt', 'w') as out:
    out.write(' '.join(str(colours.count(i)) for i in range(16)))
    out.write('\n')

