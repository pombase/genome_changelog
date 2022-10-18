#!/usr/bin/env bash
# Get the revisions where each contig was modified
contigs='chromosome1 chromosome2 chromosome3 mating_type_region pMIT'

for contig in $contigs
do
    # Create the folder structure
    mkdir -p data
    mkdir -p "data/${contig}"
    mkdir -p "data/${contig}/change_log"
    mkdir -p "data/${contig}/change_log/locations"
    mkdir -p "data/${contig}/change_log/qualifiers"

    # Download
    echo "downloading ${contig} revisions..."
    svn log -q svn+ssh://manu@curation.pombase.org/var/svn-repos/pombe-embl trunk/${contig}.contig |awk 'NR%2==0'|cut -d ' ' -f1,3,5|sed 's/^.\{1\}//' > data/${contig}/revisions.txt
done

