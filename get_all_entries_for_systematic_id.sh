echo ''
echo 'Changes in location:'
head -n 1 results/all_coordinate_changes_file_comments.tsv
grep "	$1	" results/all_coordinate_changes_file_comments.tsv
grep "	$1	" results/pre_svn_coordinate_changes_file_comments.tsv
echo ''

echo ''
echo 'Changes in qualifiers:'
head -n 1 results/all_qualifier_changes_file.tsv
grep "	$1	" results/all_qualifier_changes_file.tsv
grep "	$1	" gene_changes_comments_and_pmids/pre_svn_qualifier_changes_file.tsv
echo ''