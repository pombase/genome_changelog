# Downloads the record of genome changes that was manually kept and merges the qualifier_changes changelog (from the repository and the latest release).

set -e

mkdir -p gene_changes_comments_and_pmids/

# Coordinate changes data
echo 'systematic_id	primary_name	old_coordinates	new_coordinates	comments	reference	date' > gene_changes_comments_and_pmids/gene-coordinate-change-data.tsv
curl -k https://raw.githubusercontent.com/pombase/curation/master/data_files/gene-coordinate-change-data.tsv|tail -n+2 >> gene_changes_comments_and_pmids/gene-coordinate-change-data.tsv

# Deleted genes comments
echo 'systematic_id	comment_removal	reference_removal	date' > gene_changes_comments_and_pmids/removed-gene-data.tsv
curl -k https://raw.githubusercontent.com/pombase/curation/master/data_files/removed-gene-data.tsv|tail -n+2 >> gene_changes_comments_and_pmids/removed-gene-data.tsv

# Newly added genes comments
echo 'systematic_id	comment_addition	reference_addition	date' > gene_changes_comments_and_pmids/new-gene-data.tsv
curl -k https://raw.githubusercontent.com/pombase/curation/master/data_files/new-gene-data.tsv|tail -n+2 >> gene_changes_comments_and_pmids/new-gene-data.tsv

# qualifier_changes from pre-svn downloaded from release
curl -L https://github.com/pombase/genome_changelog/releases/download/v0.1/pre_svn_qualifier_changes_file.tsv > gene_changes_comments_and_pmids/pre_svn_qualifier_changes_file.tsv
