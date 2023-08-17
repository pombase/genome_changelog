python make_sgd_data2.py
python make_single_line_fasta.py data/*.fasta

# iterate over sorted files with a while loop

prev=""
find data -name '*_single_line.tsv' | sort | while read -r next; do

    if [ "$prev" == "" ]; then
        prev=$next
        continue
    fi

    echo $prev $next
    # # file name without extension and without single_line
    next_name=$(basename $next _single_line.tsv)
    diff --unified=3 $prev $next| grep '^[-]'|tail -n +2 |cut -c2- | awk -v h="$next_name" '{print $0 "\t" h}'> data/${next_name}_diff.tsv

    prev=$next
done

cat data/*diff.tsv|sort|uniq > all_previous_seqs.tsv
