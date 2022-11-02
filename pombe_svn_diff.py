"""
Calculate the differences between genome versions, relies on the folder structure defined in the readme.
"""
import os
from genome_functions import genome_dict_diff, build_seqfeature_dict
from formatting_functions import write_diff_to_files
import glob
from custom_biopython import SeqIO
import argparse
class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


parser = argparse.ArgumentParser(description=__doc__, formatter_class=Formatter)
parser.add_argument('--data_folder', default='data', help='folder where the analysis will be ran.')
args = parser.parse_args()
# Known errors (revision number)
skip_files = {
    'chromosome1': ['7485', '963','217'],
    'chromosome2': ['7477','1809', '1783','1395','1394','139','137','136','25','23'],
    'chromosome3': ['49'],
    }

for folder in glob.glob(f'{args.data_folder}/*'):
    output_folder = f'{folder}/change_log'
    contig_file_name = folder.replace(f'{args.data_folder}/','')
    with open(f'{folder}/revisions.txt') as f:
        revisions = f.read().splitlines()

    # Remove the known errors
    if contig_file_name in skip_files:
        revisions = [r for r in revisions if r.split()[0] not in skip_files[contig_file_name]]

    # Prepare first iteration
    old_genome_dict = None

    for i in range(len(revisions)-1):

        new_revision_list = revisions[i].split()
        old_revision_list = revisions[i+1].split()

        # Known errors in files

        if (contig_file_name , new_revision_list[0]) in skip_files:
            # Next genome will have to be read
            print(f'{folder}, known error, skipped diff {new_revision_list[0]} & {old_revision_list[0]}')
            old_genome_dict = None
            continue

        locations_output_file = f'{output_folder}/locations/{new_revision_list[0]}.tsv'
        qualifiers_output_file = f'{output_folder}/qualifiers/{new_revision_list[0]}.tsv'

        # Check if output files exist, if so skip
        if os.path.isfile(locations_output_file) and os.path.isfile(qualifiers_output_file):
            print(f'{folder}, skipped diff {new_revision_list[0]} & {old_revision_list[0]}')
            # Next genome will have to be read
            old_genome_dict = None
            continue

        print(f'{folder}, performing diff {new_revision_list[0]} & {old_revision_list[0]}')

        # Keep data from last iteration to avoid re-reading contig file
        new_genome_file = f'{folder}/{new_revision_list[0]}.contig'
        if old_genome_dict is not None:
            new_genome_dict = old_genome_dict
        else:
            # Avoids some encoding errors
            with open(new_genome_file, errors='replace') as ins:
                new_genome_dict = build_seqfeature_dict(SeqIO.read(ins,'embl'))

        old_genome_file = f'{folder}/{old_revision_list[0]}.contig'

        # Avoids some encoding errors
        with open(old_genome_file, errors='replace') as ins:
            old_genome_dict = build_seqfeature_dict(SeqIO.read(ins,'embl'))

        # Get diffs
        locations_added, locations_removed, qualifiers_added, qualifiers_removed = genome_dict_diff(new_genome_dict, old_genome_dict)
        write_diff_to_files(locations_added, locations_removed, qualifiers_added, qualifiers_removed, new_revision_list, locations_output_file, qualifiers_output_file)



