import unittest
from Bio import SeqIO
import os
import requests

from genome_functions import read_pombe_genome, make_synonym_dict

class GenomeFunctionTest(unittest.TestCase):

    def test_read_pombe_genome(self):
        """
        Test the read_pombe_genome function (mostly to ensure that the known_exception cases are applied correctly)
        """
        # Download necessary files if not downloaded already
        if not os.path.exists('test_folder'):
            os.makedirs('test_folder')

        test_file = 'test_folder/chromosome1_20090904.contig'
        if not os.path.exists(test_file):
            print('Downloading test file...')
            resp = requests.get('https://www.pombase.org/data/genome_sequence_and_features/artemis_files/OLD/20090904/chromosome1.contig')
            with open(test_file, 'wb') as out:
                out.write(resp.content)
            print('Test file downloaded')

        with open(test_file, errors='replace') as ins:
            contig_before = SeqIO.read(ins,'embl')

        # This genome has duplicated systematic_ids
        features_with_duplication = list(filter(lambda f: 'systematic_id' in f.qualifiers and len(f.qualifiers['systematic_id']) > 1, contig_before.features))
        self.assertNotEqual(len(features_with_duplication), 0)

        synonym_dictionary = make_synonym_dict('valid_ids_data/gene_IDs_names.tsv')
        contig_after = read_pombe_genome(test_file,'embl',synonym_dictionary, 'valid_ids_data/all_systematic_ids_ever.txt','valid_ids_data/known_exceptions.tsv')

        # Using this function they should be removed, but the total number of features should be the same
        features_with_duplication = list(filter(lambda f: 'systematic_id' in f.qualifiers and len(f.qualifiers['systematic_id']) > 1, contig_after.features))
        self.assertEqual(len(features_with_duplication), 0)
        self.assertEqual(len(contig_after.features), len(contig_before.features))







