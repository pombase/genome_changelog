import pandas
import argparse
import json

class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


parser = argparse.ArgumentParser(description=__doc__, formatter_class=Formatter)
parser.add_argument('--input_files', nargs='+', help='coordinate changes file (in order that you want them concatenated)')
parser.add_argument('--output_modified_coordinates',default='only_modified_coordinates.tsv', help='output tsv file')
parser.add_argument('--output_summary',default='genome_changes_summary.tsv', help='output tsv file')

args = parser.parse_args()

# Load coordinate changes
data = pandas.concat([
    pandas.read_csv(f, delimiter='\t', na_filter=False) for f in args.input_files
])

# We only consider CDS features in genes that have alleles with sequence errors
main_features_data = data[~data['feature_type'].isin(["5'UTR","3'UTR",'intron','promoter','LTR', 'misc_feature']) & data.systematic_id.str.startswith('SP')].copy()

# Make sure they are sorted properly (below, earliest_modification relies on 'added' rows being on top of 'removed' ones)
main_features_data['date'] = pandas.to_datetime(main_features_data['date'],utc=True)
main_features_data.sort_values(['date','revision','chromosome','systematic_id','feature_type','added_or_removed'],inplace=True, ascending=[False, False,True, True, True, True])
main_features_data['date'] = main_features_data['date'].dt.date

# See the columns that only differ in value and added_removed > coordinates were modified
d = main_features_data.drop(columns=['value', 'added_or_removed'])
modification_logical_index = d.duplicated(keep=False)

# See the columns that only differ in value and added_removed > coordinates were modified
modification_data = main_features_data[modification_logical_index].copy()
modification_data.to_csv(args.output_modified_coordinates, sep='\t', index=False)

## Classify genes by type of change =====================================================================================

# There were some cases in which features were removed accidentally and re-inserted some revisions
# later, so we cannot rely only on those which are added and removed in the same revision.
# Also, some are added and then edited, some are edited and then deleted.

data_subset = main_features_data[['systematic_id','chromosome', 'primary_name', 'added_or_removed']].copy()
data_subset['net_change'] = 0
data_subset['nb_changes'] = 1
data_subset.loc[data_subset['added_or_removed'] == 'added','net_change'] = 1
data_subset.loc[data_subset['added_or_removed'] == 'removed','net_change'] = -1
data_subset['earliest_modification'] = data_subset['net_change']

data_subset2 = data_subset.groupby('systematic_id', as_index=False).agg({'net_change': sum, 'nb_changes': sum, 'earliest_modification': lambda x: list(x)[-1], 'chromosome': lambda x: list(x)[0], 'primary_name': lambda x: list(x)[0]})
data_subset2['category'] = ''

data_subset2.loc[(data_subset2['net_change'] == 1) & (data_subset2['nb_changes'] == 1), 'category'] = 'added'
data_subset2.loc[(data_subset2['net_change'] == 1) & (data_subset2['nb_changes'] > 1), 'category'] = 'added_and_changed'
data_subset2.loc[(data_subset2['net_change'] == -1) & (data_subset2['nb_changes'] == 1), 'category'] = 'removed'
data_subset2.loc[(data_subset2['net_change'] == -1) & (data_subset2['nb_changes'] > 1), 'category'] = 'changed_and_removed'

data_subset2.loc[(data_subset2['net_change'] == 0) & (data_subset2['earliest_modification'] == -1), 'category'] = 'changed'
data_subset2.loc[(data_subset2['net_change'] == 0) & (data_subset2['nb_changes'] == 2) & (data_subset2['earliest_modification'] == 1), 'category'] = 'added_and_removed'
data_subset2.loc[(data_subset2['net_change'] == 0) & (data_subset2['nb_changes'] > 3) &  (data_subset2['earliest_modification'] == 1), 'category'] = 'added_changed_and_removed'

# Control for uneven changes
data_subset2['abs_change'] = abs(data_subset2['net_change'])
uneven_changes = data_subset2[data_subset2['abs_change'] > 1]
if not uneven_changes.empty:
    affected_ids = set(uneven_changes['systematic_id'])
    raise ValueError(f'The following ids have too many changes {affected_ids}')

data_subset2.drop(columns=['net_change', 'earliest_modification','nb_changes', 'abs_change'], inplace=True)

## Find merges as removal of feature + added as synonym of another  =====================================================================================
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
removed_or_merged_ids = set(data_subset2['systematic_id'][data_subset2['category'].str.contains('removed')])
removed_merged_data = main_features_data[main_features_data['systematic_id'].isin(removed_or_merged_ids)].copy()
removed_merged_data['revision'] = removed_merged_data['revision'].astype(str).copy()
synonym_data['revision'] = synonym_data['revision'].astype(str).copy()

merged_data = pandas.merge(removed_merged_data, synonym_data[['systematic_id', 'synonym_id', 'revision']], left_on=['systematic_id', 'revision'], right_on=['synonym_id', 'revision'], how='inner')
merged_data.rename(columns={'systematic_id_x': 'systematic_id', 'systematic_id_y': 'merged_into'}, inplace=True)

data_subset2 = data_subset2.merge(merged_data[['systematic_id','merged_into']], on=['systematic_id'], how='outer')

# Rename categories where merged happenned
data_subset2['category'] = data_subset2.apply(lambda row: row.category.replace('removed','merged') if not pandas.isna(row.merged_into) else row.category, axis=1)

# Add extra column indicating what types of feature the systematic_id has ever contained
extra_column_dataset = main_features_data[['systematic_id', 'feature_type']].groupby('systematic_id', as_index=False).agg({'feature_type': lambda x: ','.join(sorted(list(set(x))))})
data_subset2 = data_subset2.merge(extra_column_dataset, on='systematic_id')
data_subset2 = data_subset2[['systematic_id','chromosome','primary_name','feature_type','category','merged_into']]
data_subset2.to_csv(args.output_summary, sep='\t', index=False)
