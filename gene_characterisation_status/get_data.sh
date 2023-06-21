svn log https://curation.pombase.org/pombe-embl-repo > svn_revisions.txt

python format_revisions.py
python get_monthly_revisions.py

