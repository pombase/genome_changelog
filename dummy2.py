from genome_functions import read_pombe_genome, make_synonym_dict

synonym_dictionary = make_synonym_dict('valid_ids_data/gene_IDs_names.tsv')
contig = read_pombe_genome('pre_svn_data/chromosome1/20070531.contig','embl',synonym_dictionary, 'valid_ids_data/all_systematic_ids_ever.txt','valid_ids_data/known_exceptions.tsv')

for feature in contig.features:
    if 'systematic_id' in feature.qualifiers and feature.qualifiers['systematic_id'][0] == 'SPAC6F6.16c':
        print(feature)