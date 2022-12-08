set -e

# Download synonyms from PomBase
bash get_valid_ids.sh

# Remove possible old data
rm -rf data/*/change_log/*/*.tsv

# Download the list of revisions where the last changes were made
bash get_revisions_where_contigs_changed.sh last

# Download the contig files
python get_revisions_files.py

# Calculate differences
python pombe_svn_diff.py

# Merge them with the existing lists
python create_single_coordinate_changes_file.py --output_file temp.tsv
tail -n+2 all_coordinate_changes_file.tsv >>temp.tsv
mv temp.tsv all_coordinate_changes_file.tsv

python create_single_qualifier_changes_file.py --output_file temp.tsv
tail -n+2 all_qualifier_changes_file.tsv >>temp.tsv
mv temp.tsv all_qualifier_changes_file.tsv

# Update the changes on main features only
python get_modifications_on_main_features_only.py
