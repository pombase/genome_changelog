import unittest
import pandas
import glob
import re
from Bio import SeqIO

class PipelineTest(unittest.TestCase):

    def test_gene_summary(self):
        """
        Check that the info in genome_changes_summary.tsv is correct
        """
        # Get all systematic_id values in the last genome
        all_systematic_ids = set()
        for chromosome in glob.glob('data/*/'):
            latest_revision = sorted(glob.glob(f'{chromosome}/*.contig'), reverse=True, key= lambda x: int(re.match(r'.+\/(\d+).contig', x).groups()[0]))[0]
            this_chromosome = SeqIO.read(latest_revision,'embl')
            for feature in this_chromosome.features:
                # Only main features (mRNA seems to have been used only very few times)
                if feature.type not in ["5'UTR","3'UTR",'intron','promoter','LTR', 'misc_feature', 'mRNA','CDS_before','CDS_BEFORE','gene']:
                    if 'systematic_id' in feature.qualifiers:
                        all_systematic_ids.update(feature.qualifiers['systematic_id'])

        # load summary data
        data = pandas.read_csv('results/genome_changes_summary.tsv', sep='\t', na_filter=False)
        # Exclude the multi-CDS case
        data = data.loc[data.category != 'multi_CDS', :]
        # The ones with these categories should be present in the chromosome file
        changed_or_added = data.category.isin(['changed', 'added', 'added_and_changed'])
        # The ones that are present in the chromosome file
        present_in_contig = data.systematic_id.isin(all_systematic_ids)

        if any(changed_or_added != present_in_contig):
            msg = '\nids listed as added or changed that are not present in the contig:\n'
            msg = msg + "\n".join(data.systematic_id[changed_or_added & ~present_in_contig])
            msg = msg + '\nids listed as removed or merged that are present in the contig:\n'
            msg = msg + "\n".join(data.systematic_id[~changed_or_added & present_in_contig])
            msg = msg + "\n".join(data.category[~changed_or_added & present_in_contig])
            self.fail(msg)
        # print(len(all_systematic_ids))




