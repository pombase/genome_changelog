echo ''
echo 'Changes in location:'
head -n 1 results/all_coordinate_changes_file.tsv
grep "	$1	" results/all_coordinate_changes_file.tsv
grep "	$1	" results/pre_svn_coordinate_changes_file.tsv
echo ''

echo ''
echo 'Changes in qualifiers:'
head -n 1 results/all_qualifier_changes_file.tsv
grep "	$1	" results/all_qualifier_changes_file.tsv
grep "	$1	" results/pre_svn_qualifier_changes_file.tsv
echo ''

echo ''
echo 'Gene modifications (perhaps linked to db_xref or pombase comments):'
head -n 1 results/only_modified_coordinates_with_comments.tsv
grep "	$1	" results/only_modified_coordinates_with_comments.tsv
echo ''