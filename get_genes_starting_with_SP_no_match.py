from genome_functions import make_synonym_dict
import pandas

with open('valid_ids_data/all_systematic_ids_ever.txt') as f:
    valid_ids = set([line.strip() for line in f])
with open('valid_ids_data/all_genes_ever.txt') as f:
    all_genes = set([line.strip() for line in f])

obsoleted_ids = set(pandas.read_csv('valid_ids_data/obsoleted_ids.tsv',sep='\t')['obsolete_id'])

synonyms = make_synonym_dict('valid_ids_data/gene_IDs_names.tsv')

synonyms = set(s for s in synonyms if s.startswith('SP'))

orphan_genes = all_genes - valid_ids - synonyms - obsoleted_ids


with open('valid_ids_data/genes_starting_with_SP_no_match.txt', 'w') as out:
    for g in orphan_genes:
        out.write(g+'\n')
