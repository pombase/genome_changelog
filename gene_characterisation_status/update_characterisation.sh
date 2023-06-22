set -e

# This file is needed
mkdir -p ../valid_ids_data
echo 'systematic_id	primary_name	synonyms' > ../valid_ids_data/gene_IDs_names.tsv
curl -k https://www.pombase.org/data/names_and_identifiers/gene_IDs_names.tsv | tail -n+2 >> ../valid_ids_data/gene_IDs_names.tsv

# Delete old data
rm -rf data
rm -rf pre_svn_data

GREEN='\033[0;32m'
NC='\033[0m' # No Color

echo -e "${GREEN}downloading svn logs${NC}"
svn log https://curation.pombase.org/pombe-embl-repo > svn_revisions.txt

python format_revisions.py

echo -e "${GREEN}downloading files from svn${NC}"
python get_monthly_revisions.py

if [ -z "$(find data -mindepth 1 -print -quit)" ]; then
    echo "nothing to update, exiting"
    exit 0
fi

cd ..
echo -e "${GREEN}counting colours${NC}"
python colours_from_genome.py gene_characterisation_status/data/*

cd gene_characterisation_status

for f in data/*;
do
    revision=$(echo -n ${f}/counts.txt|cut -d'/' -f2)
    echo -ne "${revision}\t"
    cat $f/counts.txt
done > svn_counts.tsv

python merge_counts_new.py
python format_counts.py
