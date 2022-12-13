
from genome_functions import get_primary_name
from Bio.SeqIO.InsdcIO import _insdc_location_string as format_location

def format_location_change(feature, change, revision):
    return '\t'.join((
        *revision,
        feature.qualifiers['systematic_id'][0],
        get_primary_name(feature),
        feature.type,
        change,
        format_location(feature.location, None)
    ))

def format_qualifier_change(qualifier_tuple, change, revision):
    return '\t'.join((
        *revision,
        *qualifier_tuple[:4],
        change,
        qualifier_tuple[4]
    ))

def write_diff_to_files(locations_added, locations_removed, qualifiers_added, qualifiers_removed, new_revision_list, locations_output_file, qualifiers_output_file):
    # Get diffs

    # Format diffs for output
    locations_output = [format_location_change(f, 'added', new_revision_list) for f in locations_added]
    locations_output += [format_location_change(f, 'removed', new_revision_list) for f in locations_removed]
    qualifiers_output = [format_qualifier_change(q, 'added', new_revision_list) for q in qualifiers_added]
    qualifiers_output += [format_qualifier_change(q, 'removed', new_revision_list) for q in qualifiers_removed]

    # Write the output to text files

    with open(locations_output_file,'w') as out:
        if len(new_revision_list):
            out.write('\t'.join(['revision', 'user', 'date']) + '\t')
        out.write('\t'.join(['systematic_id', 'primary_name', 'feature_type', 'added_or_removed', 'value' ]) + '\n')
        out.write('\n'.join(sorted(locations_output)))
    with open(qualifiers_output_file,'w') as out:
        if len(new_revision_list):
            out.write('\t'.join(['revision', 'user', 'date']) + '\t')
        out.write('\t'.join(['systematic_id', 'primary_name', 'feature_type', 'qualifier_type', 'added_or_removed', 'value' ]) + '\n')
        out.write('\n'.join(sorted(qualifiers_output)))