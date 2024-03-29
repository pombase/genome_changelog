"""
Get a list of revisions where the genome sequence changed
"""
import os
from genome_functions import genome_sequences_are_different
import glob
import argparse
import pandas
class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass

parser = argparse.ArgumentParser(description=__doc__, formatter_class=Formatter)
parser.add_argument('--data_folders', nargs='+', default=glob.glob('data/*'), help='folders where the analysis will be ran.')
parser.add_argument('--output_file', help='output file')
args = parser.parse_args()
# Known errors (revision number)
skip_files = dict()
with open('valid_ids_data/revisions2skip.tsv') as ins:
    for line in ins:
        ls = line.strip().split('\t')
        skip_files[ls[0]] = ls[1].split(',')

all_changes = list()
for folder in args.data_folders:
    folder = os.path.normpath(folder)
    contig_file_name = os.path.basename(folder)

    with open(f'{folder}/revisions.txt') as f:
        revisions = f.read().splitlines()

    # Remove the known errors
    if contig_file_name in skip_files:
        revisions = [r for r in revisions if r.split()[0] not in skip_files[contig_file_name]]

    # Add first revision to the list if working with pre_svn data
    if 'pre_svn_data' in folder:
        first_revision_list = revisions[-1].split()
        all_changes.append([first_revision_list[0], first_revision_list[2], contig_file_name])

    for new_revision,old_revision in zip(revisions[:-1],revisions[1:]):

        new_revision_list = new_revision.split()
        old_revision_list = old_revision.split()
        new_genome_file = f'{folder}/{new_revision_list[0]}.contig'
        old_genome_file = f'{folder}/{old_revision_list[0]}.contig'

        print(f'{folder}, checking {new_revision_list[0]} & {old_revision_list[0]}')

        # Check if genome sequence has changed, in that case we use a custom SeqFeature that takes longer to load
        if genome_sequences_are_different(new_genome_file, old_genome_file):
            all_changes.append([new_revision_list[0], new_revision_list[2], contig_file_name])

output = pandas.DataFrame(all_changes, columns=['revision', 'date', 'chromosome'])
output.sort_values(['date','revision', 'chromosome'], ascending=[False, False, True],inplace=True)

def formatting_function(r):
    revision = r['revision']

    if r['date'] <= '2011-04-18':
        return f'https://www.pombase.org/data/genome_sequence_and_features/artemis_files/OLD/{revision}/{r["chromosome"]}.contig'
    else:
        if int(revision) < 374:
            return f'https://curation.pombase.org/pombe-embl-repo/{r["chromosome"]}.contig?p={revision}'
        else:
            return f'https://curation.pombase.org/pombe-embl-repo/{r["chromosome"]}.contig/trunk/?p={revision}'

output['reference'] = ''
output['link'] = ''
if not output.empty:
    output['link'] = output.apply(formatting_function, axis=1)

output.to_csv(args.output_file, sep='\t', index=False)
