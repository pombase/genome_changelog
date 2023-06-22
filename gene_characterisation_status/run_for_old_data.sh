set -e

svn log https://curation.pombase.org/pombe-embl-repo > svn_revisions.txt

python format_revisions.py
python get_monthly_revisions.py

cd ..
python colours_from_genome.py gene_characterisation_status/data/*
python get_ftp_site_files.py

cd gene_characterisation_status
python format_pre_svn_data.py

cd ..
python colours_from_genome.py gene_characterisation_status/pre_svn_data/*

# Gather the counts

for f in pre_svn_data/*;
do
    revision=$(echo -n ${f}/counts.txt|cut -d'/' -f2)
    echo -ne "${revision}\t"
    cat $f/counts.txt
done > pre_svn_counts.tsv

for f in data/*;
do
    revision=$(echo -n ${f}/counts.txt|cut -d'/' -f2)
    echo -ne "${revision}\t"
    cat $f/counts.txt
done > svn_counts.tsv

# Merge the counts
python merge_counts_old.py
python format_counts.py
