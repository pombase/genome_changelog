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
skip_files = {
    'chromosome1': ['7485', '963','217'],
    'chromosome2': ['7477','1809', '1783','1395','1394','139','137','136','25','23'],
    'chromosome3': ['49'],
    }

all_changes = list()
for folder in args.data_folders:
    folder = os.path.normpath(folder)
    contig_file_name = os.path.basename(folder)

    with open(f'{folder}/revisions.txt') as f:
        revisions = f.read().splitlines()

    # Remove the known errors
    if contig_file_name in skip_files:
        revisions = [r for r in revisions if r.split()[0] not in skip_files[contig_file_name]]


    for new_revision,old_revision in zip(revisions[:-1],revisions[1:]):

        new_revision_list = new_revision.split()
        old_revision_list = old_revision.split()
        new_genome_file = f'{folder}/{new_revision_list[0]}.contig'
        old_genome_file = f'{folder}/{old_revision_list[0]}.contig'

        print(f'{folder}, checking {new_revision_list[0]} & {old_revision_list[0]}')

        # Check if genome sequence has changed, in that case we use a custom SeqFeature that takes longer to load
        if genome_sequences_are_different(new_genome_file, old_genome_file):
            all_changes.append([old_revision_list[0], new_revision_list[0], new_revision_list[2], contig_file_name])

output = pandas.DataFrame(all_changes, columns=['previous_revision', 'new_revision', 'date', 'chromosome'])
output.sort_values(['date','new_revision', 'chromosome'], ascending=[False, False, True],inplace=True)

output.to_csv(args.output_file, sep='\t', index=False)
