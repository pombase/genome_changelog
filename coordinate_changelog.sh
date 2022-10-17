svn log -q svn+ssh://manu@curation.pombase.org/var/svn-repos/pombe-embl trunk/chromosome1.contig |awk 'NR%2==0'|cut -d ' ' -f1,3,5|cut -c2- > revisions_chromosome1.txt
