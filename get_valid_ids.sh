mkdir -p valid_ids_data
echo 'systematic_id	primary_name	synonyms' > valid_ids_data/gene_IDs_names.tsv
curl -k https://www.pombase.org/data/names_and_identifiers/gene_IDs_names.tsv | tail -n+2 >> valid_ids_data/gene_IDs_names.tsv