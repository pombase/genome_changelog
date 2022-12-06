from genome_functions import read_pombe_genome
import glob


for i,f in enumerate(glob.glob('pre_svn_data/chromosome3/*.contig')):
    print(i,f)
    try:
        read_pombe_genome(f,'embl','valid_ids_data/gene_IDs_names.tsv', 'valid_ids_data/all_systematic_ids_ever.txt','valid_ids_data/known_exceptions.tsv')
    except ValueError as err:
        print(err)
