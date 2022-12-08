echo ''
echo 'Changes in location:'
head -n 1 all_coordinate_changes_file.tsv
grep "	$1	" all_coordinate_changes_file.tsv
grep "	$1	" pre_svn_coordinate_changes_file.tsv
echo ''

echo ''
echo 'Changes in qualifiers:'
head -n 1 all_qualifier_changes_file.tsv
grep "	$1	" all_qualifier_changes_file.tsv
grep "	$1	" pre_svn_qualifier_changes_file.tsv
echo ''