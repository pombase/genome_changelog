set -e

# Download synonyms from PomBase (into folder valid_ids_data) and 
# the qualifier changelog from the pre-svn data into the same folder (this is used by get_info_from_changes.py)
bash get_data.sh

# Remove possible old data
rm -rf data/*/change_log/*/*.tsv

# Download the list of revisions where new changes were made (finds it using the latest revision values in results/all_coordinate_changes_file.tsv and results/all_qualifier_changes_file.tsv)
# they are written into data/chromosome1/revisions.txt, etc.
bash get_revisions_where_contigs_changed.sh last

# Download the contig files of those revisions with new changes
python get_revisions_files.py

# From the latest contigs of svn, load the obsolote_name field and create a tsv dictionary
python gather_obsoleted_names.py > valid_ids_data/obsoleted_ids.tsv

# Calculate differences between subsequent versions of the genome
python pombe_svn_diff.py

# Merge new diff tables with existing ones, for coordinates and qualifier changes.
python create_single_coordinate_changes_file.py --output_file temp.tsv
tail -n+2 results/all_coordinate_changes_file.tsv >>temp.tsv
mv temp.tsv results/all_coordinate_changes_file.tsv

python create_single_qualifier_changes_file.py --output_file temp.tsv
tail -n+2 results/all_qualifier_changes_file.tsv >>temp.tsv
mv temp.tsv results/all_qualifier_changes_file.tsv

# Update the file that lists changes (not additions and removals) on main features only
python get_info_from_changes.py

# Link changes in structures to changes in db_xref or pombase comments
python associate_comments_with_genome_changes.py

# List if the genome sequence has changed
python get_revisions_where_genome_sequence_changes.py --data data/* --output_file temp.tsv
tail -n+2 genome_sequence_changes.tsv >>temp.tsv
mv temp.tsv genome_sequence_changes.tsv

# Create extra files where the changes in type of RNA are not listed
python remove_rna_type_change.py
