"""
Checks that all genome files can be read in the pre_svn folder. Keeps track of the read ones not to repeat
"""

import glob
from genome_functions import read_pombe_genome, make_synonym_dict
import os


synonym_dictionary = make_synonym_dict('valid_ids_data/gene_IDs_names.tsv')

read_files = list()
if os.path.exists('test_folder/read_files.txt'):
    with open('test_folder/read_files.txt', 'r') as ins:
        read_files = map(str.strip, ins.readlines())


with open('test_folder/read_files.txt', 'a') as out:
    for file in sorted(glob.glob('pre_svn_data/*/*.contig')):
        if file not in read_files:
            print('> Reading', file)
            read_pombe_genome(file,'embl',synonym_dictionary, 'valid_ids_data/all_systematic_ids_ever.txt','valid_ids_data/known_exceptions.tsv')
            out.write(f'{file}\n')
        else:
            print('> Skipped', file)


