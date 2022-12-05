# Include synonyms

# python get_genes_starting_with_SP_no_match.py

# grep_arg=$(cat valid_ids_data/genes_starting_with_SP_no_match.txt|paste -sd "|" -)
# grep -E "/synonym|/gene" pre_svn_data/*/*.contig|grep -E "${grep_arg}"|rev|sort|rev> valid_ids_data/all_orphan_synonyms.txt

# python get_orphan_synonyms.py

# # Remove found ones from genes_starting_with_SP_no_match
grep_arg=$(cat valid_ids_data/missing_synonyms.tsv|grep -v '	$'|cut -f1 -d$'\t'|paste -sd "|" -)
grep -v -E "${grep_arg}" valid_ids_data/genes_starting_with_SP_no_match.txt >valid_ids_data/genes_starting_with_SP_no_match_temp.txt
# You can't grep into same file...
mv valid_ids_data/genes_starting_with_SP_no_match_temp.txt valid_ids_data/genes_starting_with_SP_no_match.txt
