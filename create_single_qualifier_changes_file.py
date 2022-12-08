import glob
import argparse
import pandas

class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


parser = argparse.ArgumentParser(description=__doc__, formatter_class=Formatter)
parser.add_argument('--data_folder', default='data', help='folder where the analysis will be ran.')
parser.add_argument('--output_file', help='output file (.tsv).')
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
    for f in glob.glob(f'{folder}/change_log/qualifiers/*.tsv'):
        with open(f) as ins:
            # Skip first line
            ins.readline()
            data_lines = [l.strip().split('\t') for l in ins]
            if len(data_lines):
                all_lines.extend(d[:3] + [chromosome] + d[3:] for d in data_lines)

data = pandas.DataFrame(all_lines,columns=['revision', 'user', 'date', 'chromosome', 'systematic_id', 'primary_name', 'feature_type', 'qualifier_type', 'added_or_removed', 'value'])

if not ('svn_2' in set(data['revision'])):
    data.revision = data.revision.astype(int)

data.sort_values(['revision','chromosome','systematic_id','feature_type', 'qualifier_type','added_or_removed'],inplace=True, ascending=[False,True, True, True, True, True])

data.to_csv(args.output_file,sep='\t',index=False)

