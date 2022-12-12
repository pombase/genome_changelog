import pandas
import argparse
import json

class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


parser = argparse.ArgumentParser(description=__doc__, formatter_class=Formatter)
parser.add_argument('--input_files', nargs='+', help='coordinate changes file (in order that you want them concatenated)')
parser.add_argument('--output_modified_coordinates',default='only_modified_coordinates.tsv', help='output tsv file')
parser.add_argument('--output_merged_genes',default='merged_genes.tsv', help='output tsv file')
parser.add_argument('--output_added_genes',default='added_genes.tsv', help='output tsv file')
parser.add_argument('--output_removed_genes',default='removed_genes.tsv', help='output tsv file')

args = parser.parse_args()

# Load coordinate changes
data = pandas.concat([
    pandas.read_csv(f, delimiter='\t', na_filter=False) for f in args.input_files
])

# We only consider CDS features in genes that have alleles with sequence errors
main_features_data = data[~data['feature_type'].isin(["5'UTR","3'UTR",'intron','promoter','LTR', 'misc_feature'])]

# See the columns that only differ in value and added_removed > coordinates were modified
d = main_features_data.drop(columns=['value', 'added_or_removed'])
logi = d.duplicated(keep=False)


modification_data = main_features_data[logi].copy()
modification_data.to_csv(args.output_modified_coordinates, sep='\t', index=False)

added_removed_merged_data = main_features_data[~logi].copy()

added_logi = added_removed_merged_data['added_or_removed'] == 'added'
added_data = added_removed_merged_data[added_logi]
added_data.to_csv(args.output_added_genes, sep='\t', index=False)

removed_merged_data = added_removed_merged_data[~added_logi].copy()

# Find merges as removal of feature + added as synonym of another
synonym_data = pandas.read_csv('gene_changes_comments_and_pmids/qualifier_changes.tsv',sep='\t',na_filter=False)
synonym_data = synonym_data[synonym_data['qualifier_type'] == 'synonym']

# Explode rows that contain multiple synonyms
synonym_data = synonym_data[synonym_data['qualifier_type'] == 'synonym'].copy()
def qualifier_value_to_list(value: str):
    # Convert a list like this: ('PMID:18641648', 'PMID:20118936') to a python string list using json module
    value_replaced = value.replace('(','[').replace(')',']').replace("'",'"').replace(',]',']')
    return json.loads(value_replaced)
synonym_data['value'] = synonym_data['value'].apply(qualifier_value_to_list)
synonym_data = synonym_data.explode('value')

# Exclude values that appear on added and removed (not actually changed)
d = synonym_data.drop(columns=['primary_name', 'added_or_removed'])
logi = d.duplicated(keep=False)
synonym_data = synonym_data[~logi].copy()

# Keep only the synonyms that start with SP for merge
synonym_data = synonym_data[synonym_data.value.str.startswith('SP')].copy()
synonym_data.rename(columns={'value': 'synonym_id'},inplace=True)

# Have to convert them to string to do the join for some reason
removed_merged_data['revision'] = removed_merged_data['revision'].astype(str).copy()
synonym_data['revision'] = synonym_data['revision'].astype(str).copy()

merged_data = pandas.merge(removed_merged_data, synonym_data[['systematic_id', 'synonym_id', 'revision']], left_on=['systematic_id', 'revision'], right_on=['synonym_id', 'revision'], how='inner')
merged_data.rename(columns={'systematic_id_x': 'systematic_id', 'systematic_id_y': 'merged_into'}, inplace=True)
merged_data.drop(columns=['added_or_removed', 'synonym_id'],inplace=True)

merged_data.to_csv(args.output_merged_genes, sep='\t', index=False)

removed_data = removed_merged_data[~removed_merged_data['systematic_id'].isin(merged_data['systematic_id'])].copy()

removed_data.to_csv(args.output_removed_genes, sep='\t', index=False)

