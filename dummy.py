import pandas

original_data = pandas.read_csv('all_coordinate_changes_file.tsv', sep='\t',na_filter=False)

original_data = original_data[~original_data['feature_type'].isin(["5'UTR","3'UTR",'intron','promoter','LTR', 'misc_feature'])]

sorted_data = original_data.sort_values(['systematic_id','revision','added_or_removed'],ascending=[False, False, True])
sorted_data = sorted_data.reset_index()

sorted_data.to_csv('dummy.tsv', sep='\t')
# Check if added - removed pairs come from different revisions

indexes_concerned = list()

for i in range(sorted_data.shape[0]-1):
    this_row = sorted_data.iloc[i]
    if this_row['added_or_removed'] == 'added':
        next_row = sorted_data.iloc[i+1]
        # We test equality of feature_type because sometimes that's what changed
        if (next_row['systematic_id'] == this_row['systematic_id']) and (next_row['added_or_removed'] == 'removed') and (next_row['revision'] != this_row['revision']):
            indexes_concerned.append(this_row.to_dict())
            indexes_concerned.append(next_row.to_dict())

out_data = pandas.DataFrame(indexes_concerned).drop(columns=['index'])
# out_data = out_data.sort_values(['revision','systematic_id','added_or_removed'],ascending=[False, False, True])
out_data.to_csv('dummy2.tsv',sep='\t',index=False)

# Drop those that have equal values


out_data[~out_data['value'].duplicated(keep=False)].to_csv('dummy3.tsv',sep='\t',index=False)