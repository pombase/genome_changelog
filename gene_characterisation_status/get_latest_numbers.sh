mkdir -p latest_revision

contigs='chromosome1 chromosome2 chromosome3 mating_type_region pMIT telomeric'

for contig_file in $contigs;
do
    svn export --force https://curation.pombase.org/pombe-embl-repo/trunk/${contig_file}.contig latest_revision/${contig_file}.contig
done

cd ..
echo 'systematic_id	primary_name	synonyms' > valid_ids_data/gene_IDs_names.tsv
curl -k https://www.pombase.org/data/names_and_identifiers/gene_IDs_names.tsv | tail -n+2 >> valid_ids_data/gene_IDs_names.tsv
python colours_from_genome.py gene_characterisation_status/latest_revision/
