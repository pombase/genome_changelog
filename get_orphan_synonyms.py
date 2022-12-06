import pandas
from genome_functions import read_pombe_genome


with open('valid_ids_data/all_orphan_synonyms.txt') as ins:
    all_lines = list()
    for line in ins:
        ls = line.strip().replace('\s+',' ').replace('/gene=',' gene ').replace('/synonym=',' synonym ').replace('"','').split()
        # Remove trailing FT
        ls[0] = ls[0][:-3]
        all_lines.append(ls)

data = pandas.DataFrame(all_lines,columns=['file_name','type', 'orphan_id'])

# Sort by date and type (we want the latest, and synonym if possible over gene)
data.sort_values(['orphan_id','type','file_name'], ascending=[True,False,False], inplace=True)


output_dict = dict()

for orphan_id in pandas.unique(data['orphan_id']):
    # By default, an empty string
    output_dict[orphan_id] = '\t\t'
    data_subset = data[data['orphan_id'] == orphan_id]
    for i, row in data_subset.iterrows():
        file_name = row['file_name']
        print('looking for',orphan_id,'in',file_name)

        seq_record = read_pombe_genome(file_name,'embl','valid_ids_data/gene_IDs_names.tsv', 'valid_ids_data/all_systematic_ids_ever.txt','valid_ids_data/known_exceptions.tsv')
        for feature in seq_record.features:
            if (row['type'] in feature.qualifiers) and (orphan_id in feature.qualifiers[row['type']]):
                if 'systematic_id' in feature.qualifiers:
                    output_dict[orphan_id] = ','.join(feature.qualifiers['systematic_id']) + f'\t{row["type"]}\t{file_name}'
                    break

        # If no systematic_id was found to match, try next file
        else:
            continue

        # Else move on to next orphan_id
        break

with open('valid_ids_data/missing_synonyms.tsv','w') as out:
    lines = list()
    out.write(f'orphan_id\tsystematic_id\ttype\tfound_in\n')
    for orphan_id in pandas.unique(data['orphan_id']):
        out.write(f'{orphan_id}\t{output_dict[orphan_id]}\n')
