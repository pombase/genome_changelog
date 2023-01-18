set -e

# Download synonyms from PomBase (into folder valid_ids_data)
bash get_valid_ids.sh

# Download data to link comments / changes in db_xref qualifiers with coordinate changes (into folder gene_changes_comments_and_pmids)
bash get_data_gene_changes_comments_and_pmids.sh

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
python get_info_from_changes.py --input_files all_coordinate_changes_file.tsv pre_svn_coordinate_changes_file.tsv --output_modified_coordinates only_modified_coordinates.tsv

# Link changes in structures to changes in db_xref or pombase comments
python associate_comments_with_genome_changes.py

# List if the genome has changed
python get_revisions_where_genome_sequence_changes.py --data data/* --output_file temp.tsv
tail -n+2 genome_sequence_changes.tsv >>temp.tsv
mv temp.tsv genome_sequence_changes.tsv