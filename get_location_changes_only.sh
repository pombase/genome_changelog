head -n 1 all_coordinate_changes_file.tsv
grep "	$1" all_coordinate_changes_file.tsv
grep "	$1" pre_svn_coordinate_changes_file.tsv
echo ''