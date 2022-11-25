import pandas
import json

changelog_pombase = pandas.read_csv('gene_changes_comments_and_pmids/gene-coordinate-change-data.tsv',sep='\t',na_filter=False)

changelog_script = pandas.read_csv('only_modified_coordinates.tsv',sep='\t',na_filter=False)
db_xref_script = pandas.read_csv('gene_changes_comments_and_pmids/qualifier_changes.tsv',sep='\t',na_filter=False)

changelog_script_common = changelog_script[changelog_script['systematic_id'].isin(changelog_pombase['systematic_id'])]

changelog_script_common['systematic_id'].value_counts().to_csv('counts.tsv',sep='\t')

# We find the db_xref that match a change in coordinates (same revision). To not have duplicated columns we keep the added only
matches = pandas.merge(changelog_script[changelog_script['added_or_removed'] == 'added'].drop(columns=['added_or_removed']), db_xref_script[['systematic_id', 'revision','value','added_or_removed']],on=['systematic_id', 'revision'])

# We merge multiple lines that have several additions / removals on separate lines
d = matches.drop(columns=['value_y'])
logi = d.duplicated(keep=False)

multi_row_data = matches[logi].copy()

# multi_row_data['value_fixed'] = ''

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

dbxref2include = multi_row_data.groupby(['revision','systematic_id','added_or_removed']).agg({'value_y': aggregating_function})

out = 