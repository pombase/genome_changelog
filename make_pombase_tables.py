import pandas

data = pandas.read_csv('results/genome_changes_summary_comments.tsv', sep='\t', na_filter=False)

# Remove some known errors
data = data.loc[~data.systematic_id.str.contains('controlled_curation'), :].copy()

data['comment_addition'] = data.comment_addition.apply(lambda x: '. '.join([i for i in x.split('|') if i !='']))
data['comment_removal'] = data.comment_removal.apply(lambda x: '. '.join([i for i in x.split('|') if i !='']))
data['reference_addition'] = data.reference_addition.apply(lambda x: ','.join([i for i in x.split('|') if i !='']))
data['reference_removal'] = data.reference_removal.apply(lambda x: ','.join([i for i in x.split('|') if i !='']))



added_genes = data.loc[data.category.isin(['added', 'added_and_changed']), ['systematic_id', 'primary_name', 'earliest_change', 'comment_addition', 'reference_addition']]
added_genes.rename(inplace=True,columns={
    'systematic_id': 'Systematic ID',
    'primary_name': 'Primary name',
    'earliest_change': 'Date added',
    'comment_addition': 'Comment',
    'reference_addition': 'Reference',
    })

added_genes.to_csv('results/pombase_tables/new-gene-data.tsv', sep='\t', index=False)

removed_genes = data.loc[data.category.str.contains('merged')|data.category.str.contains('removed'), ['systematic_id', 'primary_name', 'latest_change', 'latest_coords', 'comment_removal', 'reference_removal', 'merged_into']]
merged_genes = removed_genes['merged_into'] != ''

def formatting_function(r):
    if r.merged_into in r.comment_removal:
        return r.comment_removal
    elif r.comment_removal:
        return f'{r.comment_removal}. Merged into {r.merged_into}'
    else:
        return f'Merged into {r.merged_into}'

removed_genes.loc[merged_genes, 'comment_removal'] = removed_genes[merged_genes].apply(formatting_function, axis=1)
removed_genes.drop(columns='merged_into', inplace=True)
removed_genes.rename(inplace=True,columns={
    'systematic_id': 'Systematic ID',
    'primary_name': 'Primary name',
    'latest_change': 'Date removed',
    'latest_coords': 'Coordinates when removed',
    'comment_removal': 'Comment',
    'reference_removal': 'Reference',
    }
)

removed_genes.to_csv('results/pombase_tables/removed-gene-data.tsv', sep='\t', index=False)


# Changes table
data = pandas.read_csv('results/only_modified_coordinates_comments_no_type_change.tsv', sep='\t', na_filter=False)

def formatting_function(r):
    if r['db_xref']:
        if r['pombase_reference']:
            return f'{r["pombase_reference"]},{r["db_xref"]}'
        return r['db_xref']
    return r['pombase_reference']

data['reference'] = data.apply(formatting_function,axis=1)
data['location'] = data.apply(lambda r: r['chromosome'] + ':' + r['value'],axis=1)
data.rename(columns={
    'reference': 'Reference',
    'location': 'Coordinates',
    'systematic_id': 'Systematic ID',
    'primary_name' : 'Primary name',
    'added_or_removed': 'Before / After change',
    'pombase_comments': 'Comment',
    'date': 'Date'
}, inplace=True)

removed = data['Before / After change'] == 'removed'
data.loc[removed, 'Before / After change'] = 'before'
data.loc[~removed, 'Before / After change'] = 'after'

data[['Date', 'Systematic ID', 'Primary name', 'Before / After change', 'Coordinates', 'Comment', 'Reference']].to_csv('results/pombase_tables/gene-coordinate-change-data.tsv', sep='\t', index=False)
