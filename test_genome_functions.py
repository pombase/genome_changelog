import unittest
from Bio import SeqIO
import os
import requests

from genome_functions import read_pombe_genome, make_synonym_dict, get_locus_main_feature, get_locus_reference, build_seqfeature_dict

class GenomeFunctionTest(unittest.TestCase):

    def test_read_pombe_genome(self):
        """
        Test the read_pombe_genome function (mostly to ensure that the known_exception cases are applied correctly)
        """
        # Download necessary files if not downloaded already
        if not os.path.exists('test_folder'):
            os.makedirs('test_folder')

        test_file1 = 'test_folder/chromosome1_20090904.contig'
        if not os.path.exists(test_file1):
            print('Downloading test file...')
            resp = requests.get('https://www.pombase.org/data/genome_sequence_and_features/artemis_files/OLD/20090904/chromosome1.contig')
            with open(test_file1, 'wb') as out:
                out.write(resp.content)
            print('Test file 1 downloaded')

        with open(test_file1, errors='replace') as ins:
            contig_before = SeqIO.read(ins,'embl')

        # CASE 1 -> remove one of the systematic ids and keep only one
        # This genome has duplicated systematic_ids
        features_with_duplication = list(filter(lambda f: 'systematic_id' in f.qualifiers and len(f.qualifiers['systematic_id']) > 1, contig_before.features))
        self.assertNotEqual(len(features_with_duplication), 0)

        synonym_dictionary = make_synonym_dict('valid_ids_data/gene_IDs_names.tsv')
        contig_after = read_pombe_genome(test_file1,'embl',synonym_dictionary, 'valid_ids_data/all_systematic_ids_ever.txt','valid_ids_data/known_exceptions.tsv')

        # Using this function they should be removed, but the total number of features should be the same
        features_with_duplication = list(filter(lambda f: 'systematic_id' in f.qualifiers and len(f.qualifiers['systematic_id']) > 1, contig_after.features))
        self.assertEqual(len(features_with_duplication), 0)
        self.assertEqual(len(contig_after.features), len(contig_before.features))

        # CASE 2 -> Duplicate the feature (value is 'replace'). an mRNA considered to be polycistronic, so the same mRNA feature had both gene ids
        test_file2 = 'test_folder/chromosome1_20090819.contig'
        if not os.path.exists(test_file2):
            print('Downloading test file...')
            resp = requests.get('https://www.pombase.org/data/genome_sequence_and_features/artemis_files/OLD/20090819/chromosome1.contig')
            with open(test_file2, 'wb') as out:
                out.write(resp.content)
            print('Test file 2 downloaded')

        with open(test_file2, errors='replace') as ins:
            contig_before = SeqIO.read(ins,'embl')

        # This genome has a mRNA feature with these two gene qualifiers:
        features_with_duplication = list(filter(lambda f: 'gene' in f.qualifiers and ('SPAC8E11.03c' in f.qualifiers['gene'] and 'SPAC8E11.02c' in f.qualifiers['gene']), contig_before.features))
        self.assertEqual(len(features_with_duplication), 1)
        self.assertEqual(features_with_duplication[0].type, 'mRNA')

        synonym_dictionary = make_synonym_dict('valid_ids_data/gene_IDs_names.tsv')
        contig_after = read_pombe_genome(test_file2,'embl',synonym_dictionary, 'valid_ids_data/all_systematic_ids_ever.txt','valid_ids_data/known_exceptions.tsv')

        # The feature has been splitted into two
        features_without_duplication = list(filter(lambda f: f.type == 'mRNA' and 'systematic_id' in f.qualifiers and ('SPAC8E11.03c' in f.qualifiers['systematic_id'] or 'SPAC8E11.02c' in f.qualifiers['systematic_id']), contig_after.features))
        self.assertEqual(len(features_without_duplication), 2)

        # CASE 3 -> Skip the feature. This is currently not used, because in the case, the gene qualifiers never made it to be systematic_ids, but we test the implementation
        # We create a dummy all_systematic_ids_ever file with SPAC1348.14c and SPAPB8B6.01c
        custom_valid_ids_file = 'test_folder/custom_valid_ids_file.tsv'
        with open('valid_ids_data/all_systematic_ids_ever.txt') as ins:
            updated_lines = ins.readlines()
            updated_lines.append('SPAC1348.14c\n')
            updated_lines.append('SPAPB8B6.01c\n')

        with open(custom_valid_ids_file, 'w') as out:
            out.writelines(updated_lines)

        test_file3 = 'test_folder/chromosome1_20020905.contig'
        if not os.path.exists(test_file3):
            print('Downloading test file...')
            resp = requests.get('https://www.pombase.org/data/genome_sequence_and_features/artemis_files/OLD/20020905/chromosome1.contig')
            with open(test_file3, 'wb') as out:
                out.write(resp.content)
            print('Test file 3 downloaded')

        with open(test_file3, errors='replace') as ins:
            contig_before = SeqIO.read(ins,'embl')

        # This genome has a feature with these two gene qualifiers:
        features_with_duplication = list(filter(lambda f: 'gene' in f.qualifiers and ('SPAC1348.14c' in f.qualifiers['gene'] and 'SPAPB8B6.01c' in f.qualifiers['gene']), contig_before.features))
        self.assertEqual(len(features_with_duplication), 1)

        # There should be one less feature with SPAC1348.14c as gene qualifier in the before
        contig_after = read_pombe_genome(test_file3,'embl',synonym_dictionary, custom_valid_ids_file,'valid_ids_data/known_exceptions.tsv')

        features_id_before = list(filter(lambda f: f.type == 'CDS' and 'gene' in f.qualifiers and 'SPAC1348.14c' in f.qualifiers['gene'], contig_before.features))
        features_id_after = list(filter(lambda f: f.type == 'CDS' and 'gene' in f.qualifiers and 'SPAC1348.14c' in f.qualifiers['gene'], contig_after.features))

        self.assertEqual(len(features_id_before),1)
        self.assertEqual(len(features_id_after),0)

    def test_get_main_feature(self):
        test_file1 = 'test_folder/chromosome2_svn8941.contig'
        if not os.path.exists(test_file1):
            print('Downloading test file...')
            resp = requests.get('https://curation.pombase.org/pombe-embl-repo/trunk/chromosome2.contig?p=8941')
            with open(test_file1, 'wb') as out:
                out.write(resp.content)
            print('Test file downloaded')

        with open(test_file1, errors='replace') as ins:
            contig_dict = build_seqfeature_dict(SeqIO.read(ins,'embl'), False)

        # An RNA with introns
        self.assertEqual(get_locus_main_feature(contig_dict['SPBTRNALEU.06']).type, 'tRNA')

        # A CDS
        self.assertEqual(get_locus_main_feature(contig_dict['SPBC1685.17']).type, 'CDS')

    def test_get_locus_reference(self):
        test_file1 = 'test_folder/chromosome2_svn8941.contig'
        if not os.path.exists(test_file1):
            print('Downloading test file...')
            resp = requests.get('https://curation.pombase.org/pombe-embl-repo/trunk/chromosome2.contig?p=8941')
            with open(test_file1, 'wb') as out:
                out.write(resp.content)
            print('Test file downloaded')

        with open(test_file1, errors='replace') as ins:
            contig_dict = build_seqfeature_dict(SeqIO.read(ins,'embl'), False)

        # From the db_xref
        self.assertEqual(get_locus_reference(get_locus_main_feature(contig_dict['SPNCRNA.4512']))[0], 'PMID:29914874')

        # From the warning
        self.assertEqual(get_locus_reference(get_locus_main_feature(contig_dict['SPBC1685.17']))[0], 'PMID:24929437')


