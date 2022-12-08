import unittest
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord
from custom_biopython import CustomSeqFeature, SeqFeature
from genome_functions import genome_dict_diff

class CustomFeatureTest(unittest.TestCase):

    def test_custom_feature(self):

        # First a case that should make no difference
        sequence = SeqRecord('ACCGGGTTTAAAA')
        normal_feature = SeqFeature(FeatureLocation(0,5,1),type='CDS')
        custom_feature = CustomSeqFeature.from_parent(normal_feature,sequence)

        normal_feature_equal = SeqFeature(FeatureLocation(0,5,1),type='CDS')
        sequence_equal = SeqRecord('ACCGGGTTTAAAA')
        custom_feature_equal = CustomSeqFeature.from_parent(normal_feature_equal,sequence_equal)

        self.assertEqual(normal_feature, normal_feature_equal)
        self.assertEqual(custom_feature, custom_feature_equal)

        # Case where underlying sequence is different but coordinates are the same
        sequence_different = SeqRecord('AAAAAAAAAAAAAAAA')
        custom_feature_different_sequence = CustomSeqFeature.from_parent(normal_feature,sequence_different)
        self.assertNotEqual(custom_feature, custom_feature_different_sequence)

        # Case where underlying sequence is the same, but coordinates change
        normal_feature_different_coords = SeqFeature(FeatureLocation(1,6,1),type='CDS')
        custom_feature_different_coords = CustomSeqFeature.from_parent(normal_feature_different_coords,sequence_equal)
        self.assertNotEqual(custom_feature, custom_feature_different_coords)

        # Case where both underlying sequence and coordinates change, but it does not affect the feature sequence
        sequence_different2 = SeqRecord('AAAAAAAAAA')
        custom_feature_equal_though_changes = CustomSeqFeature.from_parent(normal_feature_different_coords,sequence_different2)
        self.assertNotEqual(custom_feature_different_sequence, custom_feature_equal_though_changes)
        self.assertEqual(custom_feature_different_sequence.feature_sequence, custom_feature_equal_though_changes.feature_sequence)

        # Underlying sequence is the same but coordinates change
        dict_old = {'same': {'CDS': [custom_feature]}, 'different': {'CDS': [custom_feature]}}
        dict_new = {'same': {'CDS': [custom_feature_equal]}, 'different': {'CDS': [custom_feature_different_coords]}}
        locations_added, locations_removed, qualifiers_added, qualifiers_removed = genome_dict_diff(dict_new, dict_old)
        self.assertEqual(locations_added[0],custom_feature_different_coords)
        self.assertEqual(locations_removed[0],custom_feature)

        # Underlying sequence is different but coordinates are the same
        dict_old = {'same': {'CDS': [custom_feature]}, 'different': {'CDS': [custom_feature]}}
        dict_new = {'same': {'CDS': [custom_feature_equal]}, 'different': {'CDS': [custom_feature_different_sequence]}}
        locations_added, locations_removed, qualifiers_added, qualifiers_removed = genome_dict_diff(dict_new, dict_old)
        self.assertEqual(locations_added[0],custom_feature_different_sequence)
        self.assertEqual(locations_removed[0],custom_feature)

        # Sequence and coordinates change, but feature sequence doesn't
        dict_old = {'same': {'CDS': [custom_feature]}, 'different': {'CDS': [custom_feature_different_sequence]}}
        dict_new = {'same': {'CDS': [custom_feature_equal]}, 'different': {'CDS': [custom_feature_equal_though_changes]}}
        locations_added, locations_removed, qualifiers_added, qualifiers_removed = genome_dict_diff(dict_new, dict_old)
        self.assertEqual(locations_added,[])
        self.assertEqual(locations_removed,[])
