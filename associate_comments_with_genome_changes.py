import pandas
import json

# Load data ===================================================

# Data from pombase comments
changelog_pombase = pandas.read_csv('gene_changes_comments_and_pmids/gene-coordinate-change-data.tsv',sep='\t',na_filter=False)
# Fill in empty dates with the closest
changelog_pombase['date'][changelog_pombase['date'] == ''] = pandas.NaT
changelog_pombase['date'] = changelog_pombase['date'].fillna(method='bfill')

# Data generated from get_modifications_on_main_features_only.py
changelog_script = pandas.read_csv('only_modified_coordinates.tsv',sep='\t',na_filter=False)
changelog_script['original_index'] = changelog_script.index

# Data generated from 
db_xref_script = pandas.read_csv('gene_changes_comments_and_pmids/qualifier_changes.tsv',sep='\t',na_filter=False)
db_xref_script.rename(columns={'value': 'db_xref'}, inplace=True)
# We find the db_xref that match a change in coordinates (same revision). To not have duplicated columns we keep the added only
matches = pandas.merge(changelog_script[changelog_script['added_or_removed'] == 'added'].drop(columns=['added_or_removed']), db_xref_script[['systematic_id', 'revision', 'feature_type', 'db_xref','added_or_removed']],on=['systematic_id', 'revision', 'feature_type'])

matches.to_csv('matches.tsv',sep='\t', index=False)

# We merge multiple lines that have several additions / removals on separate lines
d = matches.drop(columns=['db_xref'])
logi = d.duplicated(keep=False)

multi_row_data = matches[logi].copy()

unique_identifier_cols = ['systematic_id', 'revision','added_or_removed','feature_type']

def aggregating_function(values: str):
    out_set = set()
    for value in values:
        # Value would be a string in a list in this format ('PMID:18641648', 'PMID:20118936'), we can read it using the json module
        value_replaced = value.replace('(','[').replace(')',']').replace("'",'"').replace(',]',']')
        out_set.update(json.loads(value_replaced))

    # Reformat in the same format again
    # Add the single quotes
    out_set = [f'\'{i}\'' for i in out_set]

    return f'({",".join(out_set)})'

dbxref2include = multi_row_data.groupby(unique_identifier_cols).agg({'db_xref': aggregating_function})

data2merge = matches.merge(dbxref2include,on=unique_identifier_cols,how='outer')

replace_logical_index =  ~pandas.isna(data2merge['db_xref_y'])

data2merge.loc[replace_logical_index,'db_xref_x'] = data2merge['db_xref_y'][replace_logical_index]
data2merge.drop(columns=['db_xref_y'],inplace=True)
data2merge.rename(columns={'db_xref_x': 'db_xref'}, inplace=True)

output_data = changelog_script.merge(data2merge[unique_identifier_cols + ['db_xref']], on=unique_identifier_cols,how='outer')

output_data['date'] = pandas.to_datetime(output_data['date'],utc=True)
changelog_pombase['date'] = pandas.to_datetime(changelog_pombase['date'],utc=True)

output_data = output_data.sort_values(['date','systematic_id', 'added_or_removed'],ascending=[True, False, True])
changelog_pombase = changelog_pombase.sort_values(['date'],ascending=[True])

# A note on this: merge_asof(left, right) finds for every row in left the nearest row in right, so to have a single match between
# comments in changelog_pombase to output_data rows, we have to do it like this
temp_data = pandas.merge_asof(changelog_pombase[['systematic_id','reference','comments', 'date']],output_data, by=['systematic_id'], on=['date'], direction='nearest')
temp_data.to_csv('random.tsv',sep='\t', index=False)
output_data = pandas.merge(output_data,temp_data[['original_index', 'reference','comments']], on='original_index',how='outer')

output_data['date'] = output_data['date'].dt.date

output_data = output_data.sort_values(['date','revision','systematic_id', 'added_or_removed'],ascending=[False, False, True, True])

output_data.rename(columns={'reference': 'pombase_reference'}, inplace=True)
output_data.rename(columns={'comments': 'pombase_comments'}, inplace=True)

output_data = output_data.drop_duplicates()
output_data.sort_values(['original_index'], inplace=True)
output_data.drop(columns=['original_index']).to_csv('new_matches.tsv',sep='\t', index=False)


