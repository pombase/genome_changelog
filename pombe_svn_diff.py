"""
Calculate the differences between genome versions, relies on the folder structure defined in the readme.
"""
import os
from genome_functions import genome_dict_diff, build_seqfeature_dict,read_pombe_genome, genome_sequences_are_different, make_synonym_dict
from formatting_functions import write_diff_to_files
import glob
import argparse
class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


parser = argparse.ArgumentParser(description=__doc__, formatter_class=Formatter)
parser.add_argument('--data_folders', nargs='+', default=glob.glob('data/*'), help='folders where the analysis will be ran.')
args = parser.parse_args()
# Known errors (revision number)
skip_files = {
    'chromosome1': ['8645','7485', '963','217'],
    'chromosome2': ['8648','8645','7477','1809', '1783','1395','1394','139','137','136','25','23'],
    'chromosome3': ['8645','49'],
    }

synonym_dict = make_synonym_dict('valid_ids_data/gene_IDs_names.tsv', 'valid_ids_data/obsoleted_ids.tsv', 'valid_ids_data/missing_synonyms.tsv')

for folder in args.data_folders:
    folder = os.path.normpath(folder)
    output_folder = f'{folder}/change_log'
    contig_file_name = os.path.basename(folder)

    with open(f'{folder}/revisions.txt') as f:
        revisions = f.read().splitlines()

    # Remove the known errors
    if contig_file_name in skip_files:
        revisions = [r for r in revisions if r.split()[0] not in skip_files[contig_file_name]]

    # Prepare first iteration
    old_genome_dict = None
    sequences_different_prev_iteration = None

    for i in range(len(revisions)-1):

        new_revision_list = revisions[i].split()
        old_revision_list = revisions[i+1].split()
        new_genome_file = f'{folder}/{new_revision_list[0]}.contig'
        old_genome_file = f'{folder}/{old_revision_list[0]}.contig'

        locations_output_file = f'{output_folder}/locations/{new_revision_list[0]}.tsv'
        qualifiers_output_file = f'{output_folder}/qualifiers/{new_revision_list[0]}.tsv'

        # Check if output files exist, if so skip
        if os.path.isfile(locations_output_file) and os.path.isfile(qualifiers_output_file):
            print(f'{folder}, skipped diff {new_revision_list[0]} & {old_revision_list[0]}')
            # Next genome will have to be read
            old_genome_dict = None
            continue

        print(f'{folder}, performing diff {new_revision_list[0]} & {old_revision_list[0]}')

        # Check if genome sequence has changed, in that case we use a custom SeqFeature that takes longer to load
        sequences_different_this_iteration = genome_sequences_are_different(new_genome_file, old_genome_file)

        if (sequences_different_this_iteration):
            print('> Loading feature sequences, this comparison will be slow')

        # Keep data from last iteration to avoid re-reading contig file
        if (old_genome_dict is not None) and (sequences_different_this_iteration == sequences_different_prev_iteration):
            new_genome_dict = old_genome_dict
        else:
            new_genome_dict = build_seqfeature_dict(read_pombe_genome(new_genome_file, 'embl', synonym_dict, 'valid_ids_data/all_systematic_ids_ever.txt','valid_ids_data/known_exceptions.tsv'), sequences_different_this_iteration)

        old_genome_dict = build_seqfeature_dict(read_pombe_genome(old_genome_file, 'embl', synonym_dict, 'valid_ids_data/all_systematic_ids_ever.txt','valid_ids_data/known_exceptions.tsv'), sequences_different_this_iteration)

        # Get diffs
        locations_added, locations_removed, qualifiers_added, qualifiers_removed = genome_dict_diff(new_genome_dict, old_genome_dict,sequences_different_this_iteration)
        write_diff_to_files(locations_added, locations_removed, qualifiers_added, qualifiers_removed, new_revision_list, locations_output_file, qualifiers_output_file)



