"""
Compare two genomes, vanilla version (relies exclusively on systematic_id + feature_type as a unique identifier for a feature
and does not use synonyms nor deals witch changes in underlying genome sequence). To be used carefully.
"""

import argparse
from genome_functions import genome_dict_diff, build_seqfeature_dict
from formatting_functions import write_diff_to_files
from custom_biopython import SeqIO

class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass

parser = argparse.ArgumentParser(description=__doc__, formatter_class=Formatter)
parser.add_argument('--new_genome', type=str)
parser.add_argument('--old_genome', type=str)
parser.add_argument('--revision_string', type=str, default='', help='Revision identifier associated with the diff. In the example cases, the svn revision or FTP folder name.')
parser.add_argument('--format', type=str, default='embl', help='Format of the genomes (to be passed to SeqIO.read)')
parser.add_argument('--output_locations_file', type=str)
parser.add_argument('--output_qualifiers_file', type=str)

args = parser.parse_args()

with open(args.new_genome, errors='replace') as ins:
    new_genome_dict = build_seqfeature_dict(SeqIO.read(ins,args.format), False)

with open(args.old_genome, errors='replace') as ins:
    old_genome_dict = build_seqfeature_dict(SeqIO.read(ins,args.format), False)

locations_added, locations_removed, qualifiers_added, qualifiers_removed = genome_dict_diff(new_genome_dict, old_genome_dict, False)
write_diff_to_files(locations_added, locations_removed, qualifiers_added, qualifiers_removed, args.revision_string.split(), args.output_locations_file, args.output_qualifiers_file)