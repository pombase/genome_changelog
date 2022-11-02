contigs='chromosome1 chromosome2 chromosome3 mating_type_region pMIT'

for contig in $contigs
do
    # Create the folder structure
    mkdir -p pre_svn_data
    mkdir -p "pre_svn_data/${contig}"
    mkdir -p "pre_svn_data/${contig}/change_log"
    mkdir -p "pre_svn_data/${contig}/change_log/locations"
    mkdir -p "pre_svn_data/${contig}/change_log/qualifiers"

    # Download files from first revision of svn
    svn cat -r 2 svn+ssh://manu@curation.pombase.org/var/svn-repos/pombe-embl/trunk/${contig}.contig > "pre_svn_data/${contig}/svn_2.contig"
done