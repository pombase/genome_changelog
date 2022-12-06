# Stop on errors
set -e

# From the pre_svn contigs get all the /systematic_id qualifiers
grep --no-filename '/systematic_id' pre_svn_data/*/*.contig | sort | uniq|cut -c38-|rev | cut -c 2- | rev > valid_ids_data/all_systematic_ids_ever.txt

# From the pre_svn contigs get all the /gene qualifiers in which the value starts by SP
grep --no-filename '/gene' pre_svn_data/*/*.contig | sort | uniq|cut -c29-|rev | cut -c 2- | rev|grep SP|grep -v '=' > valid_ids_data/all_genes_starting_with_SP_ever.txt

# From the latest contigs of svn, load the obsolote_name field and create a tsv dictionary
python gather_obsoleted_names.py > valid_ids_data/obsoleted_ids.tsv

# Get all identifiers that start with SP but do not match any current id or synonym and write them to genes_starting_with_SP_no_match.txt
python get_genes_starting_with_SP_no_match.py

# Find those identifiers in \synonym or \gene qualifiers. If any of those qualifiers is in a feature that has a \gene or \systematic_id corresponding to a current
# systematic id, assign the unknown identifier as a synonym. Prioritises the use of /synonym qualifiers over /gene qualifiers.
grep_arg=$(cat valid_ids_data/genes_starting_with_SP_no_match.txt|awk '{print "\"" $0 "\""}'|paste -sd '|' -)
grep -E "/synonym|/gene" pre_svn_data/*/*.contig|grep -E "${grep_arg}"|rev|sort|rev> valid_ids_data/all_orphan_synonyms.txt
python get_orphan_synonyms.py

# Remove found ones from genes_starting_with_SP_no_match
grep_arg=$(cat valid_ids_data/missing_synonyms.tsv|grep -v '	$'|cut -f1 -d$'\t'|awk '{print "^" $0 "$"}'|paste -sd "|" -)
grep -v -E "${grep_arg}" valid_ids_data/genes_starting_with_SP_no_match.txt >valid_ids_data/genes_starting_with_SP_no_match_temp.txt
# You can't grep into same file...
mv valid_ids_data/genes_starting_with_SP_no_match_temp.txt valid_ids_data/genes_starting_with_SP_no_match.txt

# Now genes_starting_with_SP_no_match.txt contains the ones that don't match anything.

# Remove missing synonyms, which are committed to the repo
cat valid_ids_data/missing_synonyms.tsv|grep -v '	$' > valid_ids_data/missing_synonyms_temp.tsv
mv valid_ids_data/missing_synonyms_temp.tsv valid_ids_data/missing_synonyms.tsv