"""
Create a new file where only unique remove coordinates are shown. Comments and db_xref are combined comma separated.
PMID:xxxxxxxx values are removed
"""

import pandas

data = pandas.concat([
    pandas.read_csv(f, delimiter='\t', na_filter=False, dtype=str) for f in ['results/all_coordinate_changes_file_comments_no_type_change.tsv','results/pre_svn_coordinate_changes_file_comments_no_type_change.tsv']
])

def agg_function (x):
    """Used to aggregate rows as comma-separatated values"""
    return ','.join([ele for ele in x if ele != ''])


# Sometimes PMID:xxxxxxxx has been used, we remove those values
data.db_xref.loc[data.db_xref.str.contains('PMID:xx')] = data.db_xref.loc[data.db_xref.str.contains('PMID:xx')].apply(lambda x: ','.join( ele for ele in x.split(',') if not ele.startswith('PMID:x')))

# We also want to keep the db_xref value of added changes, we keep a separate dataset to merge
added_rows = data.loc[data.added_or_removed == 'added', ['systematic_id', 'value', 'db_xref', 'pombase_reference', 'pombase_comments']].copy()
# Aggregate the db_xref col
added_rows = added_rows.groupby(['systematic_id', 'value'], as_index=False).agg({'db_xref': agg_function, 'pombase_reference': agg_function, 'pombase_comments': agg_function})

# We keep only the removed rows (should contain all non-current values)
data = data[data.added_or_removed == 'removed'].copy()
data.drop(columns=['added_or_removed', 'revision', 'user'], inplace=True)

# Aggregate db_xref and comments, and keep only latest date
agg_cols = ['db_xref',	'pombase_reference', 'pombase_comments', 'date']
grouping_cols = list(data.columns)
[grouping_cols.remove(col) for col in agg_cols]
data = data.groupby(grouping_cols, as_index=False).agg({'date': lambda x: list(x)[-1], 'db_xref': agg_function, 'pombase_reference': agg_function, 'pombase_comments': agg_function })

# Merge db_xref of added rows
data.rename(columns={'db_xref': 'db_xref_removed', 'pombase_reference': 'pombase_reference_removed', 'pombase_comments': 'pombase_comments_removed'}, inplace=True)
added_rows.rename(columns={'db_xref': 'db_xref_added', 'pombase_reference': 'pombase_reference_added', 'pombase_comments': 'pombase_comments_added'}, inplace=True)
data = data.merge(added_rows, on=['systematic_id', 'value'], how='left')

# Re-sort columns and save
data = data[['chromosome','systematic_id','primary_name','feature_type','value','date','db_xref_removed','db_xref_added','pombase_reference_removed','pombase_reference_added','pombase_comments_removed','pombase_comments_added']]
data.to_csv('results/all_previous_coords.tsv', sep='\t', index=False)
