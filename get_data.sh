mkdir -p valid_ids_data
echo 'systematic_id	primary_name	synonyms' > valid_ids_data/gene_IDs_names.tsv
curl -k https://www.pombase.org/data/names_and_identifiers/gene_IDs_names.tsv | tail -n+2 >> valid_ids_data/gene_IDs_names.tsv

# qualifier_changes from pre-svn downloaded from release
curl -L https://github.com/pombase/genome_changelog/releases/latest/download/pre_svn_qualifier_changes_file.tsv > gene_changes_comments_and_pmids/pre_svn_qualifier_changes_file.tsv