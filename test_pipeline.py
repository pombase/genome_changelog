import unittest
from tempfile import TemporaryDirectory
from associate_comments_with_genome_changes import main as associate_comments
import pandas
class PipelineTest(unittest.TestCase):

    def test_associate_comments(self):
        """
        Ensure that the script associate_comments_with_genome_changes does not alter any value in pre-existing columns, or change the number of rows.
        """

        with TemporaryDirectory() as d:
            associate_comments(f'{d}/comments_associated.tsv')
            new_data = pandas.read_csv(f'{d}/comments_associated.tsv', sep ='\t', na_filter=False)
            old_data = pandas.read_csv(f'only_modified_coordinates.tsv', sep ='\t', na_filter=False)
            new_data.drop(columns=['db_xref', 'pombase_reference', 'pombase_comments'], inplace=True)
            self.assertEqual(new_data.shape,old_data.shape,'The shapes of the datasets differ. Likely because associate_comments_with_genome_changes adds additional rows')


