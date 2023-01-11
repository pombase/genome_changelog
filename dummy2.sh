# extra CDS at rev 196

for id in $(cat problem_list.tsv)
do
    bash get_location_changes_only.sh $id
done