import glob
import argparse
class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


parser = argparse.ArgumentParser(description=__doc__, formatter_class=Formatter)
parser.add_argument('--data_folder', default='data', help='folder where the analysis will be ran.')
args = parser.parse_args()

chromosome_dictionary = {
    'chromosome1': 'I',
    'chromosome2': 'II',
    'chromosome3': 'III',
    'mating_type_region': 'mating_type_region',
    'pMIT': 'mitochondrial',
    }

all_lines = list()
for folder in glob.glob(f'{args.data_folder}/*'):
    contig_file_name = folder.replace(f'{args.data_folder}','').replace('/','')
    chromosome = chromosome_dictionary[contig_file_name]
    for f in glob.glob(f'{folder}/change_log/locations/*.tsv'):
        with open(f) as ins:
            # Skip first line
            ins.readline()
            data_lines = [l.strip().split('\t') for l in ins]
            if len(data_lines):
                d = data_lines[0]
            all_lines.extend('\t'.join(d[:3] + [chromosome] + d[3:]) for d in data_lines)

print('revision', 'user', 'date', 'chromosome', 'systematic_id', 'primary_name', 'feature_type', 'added_or_removed', 'value',sep='\t')
if len(all_lines):
    if args.data_folder == 'data':
        sorting_fun = lambda l: int(l.split('\t')[0])
    else:
        sorting_fun = lambda l: l
    print(*sorted(all_lines, reverse=True, key=sorting_fun),sep='\n')



