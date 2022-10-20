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
    mkdir -p "data/${contig}/diff"

    last_revision_coordinates=$(sed -n '2p' all_coordinate_changes_file.tsv|cut -d$'\t' -f1)
    last_revision_qualifiers=$(sed -n '2p' all_qualifier_changes_file.tsv|cut -d$'\t' -f1)

    # Get the maximum of the two
    last_revision_included=$(( $last_revision_coordinates > $last_revision_qualifiers ? $last_revision_coordinates : $last_revision_qualifiers ))

    # Download
    echo "downloading ${contig} revisions from ${last_revision_included}..."
    if [ $1 == 'all' ];then
        svn log -q svn+ssh://manu@curation.pombase.org/var/svn-repos/pombe-embl trunk/${contig}.contig |awk 'NR%2==0'|cut -d ' ' -f1,3,5|sed 's/^.\{1\}//' > data/${contig}/revisions.txt
    elif [ $1 == 'last' ];then
        svn log -q -r HEAD:${last_revision_included}  svn+ssh://manu@curation.pombase.org/var/svn-repos/pombe-embl trunk/${contig}.contig |awk 'NR%2==0'|cut -d ' ' -f1,3,5|sed 's/^.\{1\}//' > data/${contig}/revisions.txt
    fi

done

