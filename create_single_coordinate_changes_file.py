import glob

skip_files = {
    'chromosome1': 'I',
    'chromosome2': 'II',
    'chromosome3': 'III',
    'mating_type_region': 'mating_type_region',
    'pMIT': 'mitochondrial',
    }

all_lines = list()
for folder in glob.glob('data/*'):
    contig_file_name = folder.replace('data/','')
    chromosome = skip_files[contig_file_name]
    for f in glob.glob(f'{folder}/change_log/locations/*.tsv'):
        with open(f) as ins:
            # Skip first line
            ins.readline()
            data_lines = [l.strip().split('\t') for l in ins]
            if len(data_lines):
                d = data_lines[0]
            all_lines.extend('\t'.join(d[:3] + [chromosome] + d[3:]) for d in data_lines)

print('revision', 'user', 'date', 'chromosome', 'systematic_id', 'primary_name', 'feature_type', 'added_or_removed', 'value',sep='\t')
print(*sorted(all_lines, reverse=True, key= lambda l: int(l.split('\t')[0])),sep='\n')



