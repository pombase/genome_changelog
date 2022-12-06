"""
Downloads the record of genome changes that was manually kept and merges the qualifier_changes changelog (from the repository and the latest release).
"""

mkdir -p gene_changes_comments_and_pmids/

# rename the columns
echo 'systematic_id	primary_name	old_coordinates	new_coordinates	comments	reference	date' > gene_changes_comments_and_pmids/gene-coordinate-change-data.tsv
curl -k https://raw.githubusercontent.com/pombase/curation/master/data_files/gene-coordinate-change-data.tsv|tail -n+2 >> gene_changes_comments_and_pmids/gene-coordinate-change-data.tsv

echo 'revision	user	date	chromosome	systematic_id	primary_name	feature_type	qualifier_type	added_or_removed	value' > gene_changes_comments_and_pmids/qualifier_changes.tsv

# qualifier_changes from svn downloaded from release
curl -L https://github.com/pombase/genome_changelog/releases/download/v0.1/all_qualifier_changes_file.tsv|grep '	db_xref	' >> gene_changes_comments_and_pmids/qualifier_changes.tsv

# qualifier_changes from pre-svn downloaded from release
curl -L https://github.com/pombase/genome_changelog/releases/download/v0.1/pre_svn_qualifier_changes_file.tsv|grep '	db_xref	' >> gene_changes_comments_and_pmids/qualifier_changes.tsv

# merge with changes that have not included in a release yet
cd gene_changes_comments_and_pmids/
python merge_qualifier_changes.py
