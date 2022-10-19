#%%
import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
import glob
from Bio.SeqIO.InsdcIO import _insdc_location_string as format_location
from Bio.GenBank.Scanner import EmblScanner
import re

# We override this method to allow no space between number and BP
@staticmethod
def permissive_seq_length_scanner(consumer, text):
    length_parts = text.split()
    if len(length_parts) == 1:
        length_parts = re.match('(\d+)([^\d]+)',text).groups()
    assert len(length_parts) == 2, "Invalid sequence length string %r" % text
    assert length_parts[1].upper() in ["BP", "BP.", "AA", "AA."]
    consumer.size(length_parts[0])

EmblScanner._feed_seq_length = permissive_seq_length_scanner


def features_are_equal(self, other):
    if not isinstance(other, SeqFeature):
        return False
    return (
        self.location==other.location and
        self.type==other.type and
        self.location_operator==other.location_operator and
        self.strand==other.strand and
        self.id==other.id and
        self.qualifiers==other.qualifiers and
        self.ref==other.ref and
        self.ref_db==other.ref_db
    )

## Override equality test
SeqFeature.__eq__ = features_are_equal


def build_seqfeature_dict(genome: SeqRecord):
    """
    Creates a "feature dictionary" where they keys are the systematic_id, and the values are
    "gene dictionaries". In the "gene dictionary" the keys are the types of features
    (e.g. CDS), and the values are lists containing all occurrences of that type of
    feature in the gene. E.g. in pombe:

    >> feature_dictionary = build_seqfeature_dict('chromosome1.contig')
    >> feature_dictionary['SPAC22E12.16c']['CDS'] # A list with 1 element containing the CDS feature
    >> feature_dictionary['SPAC22E12.16c']['intron'] # A list with 3 elements containing the intron features
    """
    out_dict = dict()

    for feature in genome.features:
        feature: SeqFeature
        if 'systematic_id' not in feature.qualifiers:
            continue
        gene_id = feature.qualifiers['systematic_id'][0]
        if gene_id not in out_dict:
            out_dict[gene_id] = dict()
        feature_type = feature.type
        if feature_type not in out_dict[gene_id]:
            out_dict[gene_id][feature_type] = list()
        out_dict[gene_id][feature_type].append(feature)


    return out_dict

def get_primary_name(f: SeqFeature):
    return f.qualifiers['primary_name'][0] if 'primary_name' in f.qualifiers else ''

def genome_dict_diff(new_genome_dict, old_genome_dict) -> tuple[list[SeqFeature],list[SeqFeature],list,list]:
    """
    Takes two genome dictionaries as input, returns:
    - locations_added: list of SeqFeatures in new_genome_dict that have been added or their location has changed.
    - locations_removed: list of SeqFeatures in old_genome_dict that have been removed or their location has changed.
    - qualifiers_added: list of qualifiers that have been added or modified for (systematic_id, feature_type) pairs that existed in old_genome_dict.
    - qualifiers_removed: list of qualifiers that have been removed or modified for (systematic_id, feature_type) pairs that remain in new_genome_dict.
    The qualifiers are stored as tuples of 5 elements that contain the following elements:
        - systematic_id, e.g. SPAC227.08c
        - primary_name (if exists else empty string)
        - feature_type, e.g. CDS
        - qualifier_key, e.g. product
        - qualifier_values, e.g. ('mRNA cleavage and polyadenylation specificity factor complex zinc finger subunit Yth1',)
    """
    locations_added = list()
    locations_removed = list()
    qualifiers_added = list()
    qualifiers_removed = list()

    ## First - genes that have been entirely added or removed
    # Only in the new
    for systematic_id in set(new_genome_dict.keys()) - set(old_genome_dict.keys()):
        locations_added += sum(new_genome_dict[systematic_id].values(),[])

    # Only in the old
    for systematic_id in set(old_genome_dict.keys()) - set(new_genome_dict.keys()):
        locations_removed += sum(old_genome_dict[systematic_id].values(),[])

    ## Second - existing genes that have been modified
    for systematic_id in set(old_genome_dict.keys()).intersection(set(new_genome_dict.keys())):
        new_annotation = new_genome_dict[systematic_id]
        old_annotation = old_genome_dict[systematic_id]
        if old_annotation != new_annotation:

            # A new feature_type has been added
            for feature_type in set(new_annotation.keys()) - set(old_annotation.keys()):
                locations_added += new_annotation[feature_type]

            # A feature_type has been removed
            for feature_type in set(old_annotation.keys()) - set(new_annotation.keys()):
                locations_removed += old_annotation[feature_type]

            # A feature that is in both the new and old has changed
            for feature_type in set(old_annotation.keys()).intersection(new_annotation.keys()):

                # The list of features of that feature_type for that systematic_id
                old_features = old_annotation[feature_type]
                new_features = new_annotation[feature_type]

                ## First - changes to location
                old_locations = [f.location for f in old_features]
                new_locations = [f.location for f in new_features]

                locations_added += [f for f in new_features if f.location not in old_locations]
                locations_removed += [f for f in old_features if f.location not in new_locations]

                ## Second - changes to qualifiers
                old_qualifiers = set()
                new_qualifiers = set()

                for f in old_features:
                    old_qualifiers.update( set([(systematic_id, get_primary_name(f), feature_type, key, tuple(value)) for key, value in f.qualifiers.items()]))
                for f in new_features:
                    new_qualifiers.update( set([(systematic_id, get_primary_name(f), feature_type, key, tuple(value)) for key, value in f.qualifiers.items()]))

                qualifiers_added.extend(new_qualifiers - old_qualifiers)
                qualifiers_removed.extend(old_qualifiers - new_qualifiers)

    return locations_added, locations_removed, qualifiers_added, qualifiers_removed


def format_location_change(feature: SeqFeature, change, revision):
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
        str(qualifier_tuple[4])
    ))

# Known errors
skip_files = {
    'chromosome1': ['7485', '963','217'],
    'chromosome2': ['7477','1809', '1783','1395','1394','139','137','136','25','23'],
    'chromosome3': ['49'],
    }

for folder in glob.glob('data/*'):
    output_folder = f'{folder}/change_log'
    contig_file_name = folder.replace('data/','')
    with open(f'{folder}/revisions.txt') as f:
        revisions = f.read().splitlines()

    # Remove the known errors
    if contig_file_name in skip_files:
        revisions = [r for r in revisions if r.split()[0] not in skip_files[contig_file_name]]

    # Prepare first iteration
    old_genome_dict = None

    for i in range(len(revisions[1:])-1):

        new_revision_list = revisions[i].split()
        old_revision_list = revisions[i+1].split()

        # Known errors in files

        if (contig_file_name , new_revision_list[0]) in skip_files:
            # Next genome will have to be read
            print(f'{folder}, known error, skipped diff {new_revision_list[0]} & {old_revision_list[0]}')
            old_genome_dict = None
            continue

        locations_output_file = f'{output_folder}/locations/{new_revision_list[0]}.tsv'
        qualifiers_output_file = f'{output_folder}/qualifiers/{new_revision_list[0]}.tsv'

        # Check if output files exist, if so skip
        if os.path.isfile(locations_output_file) and os.path.isfile(qualifiers_output_file):
            print(f'{folder}, skipped diff {new_revision_list[0]} & {old_revision_list[0]}')
            # Next genome will have to be read
            old_genome_dict = None
            continue

        print(f'{folder}, performing diff {new_revision_list[0]} & {old_revision_list[0]}')

        # Keep data from last iteration to avoid re-reading contig file
        new_genome_file = f'{folder}/{new_revision_list[0]}.contig'
        if old_genome_dict is not None:
            new_genome_dict = old_genome_dict
        else:
            # Avoids some encoding errors
            with open(new_genome_file, errors='replace') as ins:
                new_genome_dict = build_seqfeature_dict(SeqIO.read(ins,'embl'))

        old_genome_file = f'{folder}/{old_revision_list[0]}.contig'

        # Avoids some encoding errors
        with open(old_genome_file, errors='replace') as ins:
            old_genome_dict = build_seqfeature_dict(SeqIO.read(ins,'embl'))

        # Get diffs
        locations_added, locations_removed, qualifiers_added, qualifiers_removed = genome_dict_diff(new_genome_dict, old_genome_dict)

        # Format diffs for output
        locations_output = [format_location_change(f, 'added', new_revision_list) for f in locations_added]
        locations_output += [format_location_change(f, 'removed', new_revision_list) for f in locations_removed]
        qualifiers_output = [format_qualifier_change(q, 'added', new_revision_list) for q in qualifiers_added]
        qualifiers_output += [format_qualifier_change(q, 'removed', new_revision_list) for q in qualifiers_removed]

        # Write the output to text files

        with open(locations_output_file,'w') as out:
            out.write('\t'.join(['revision', 'user', 'date', 'systematic_id', 'primary_name', 'feature_type', 'added_or_removed', 'value' ]) + '\n')
            out.write('\n'.join(sorted(locations_output)))
        with open(qualifiers_output_file,'w') as out:
            out.write('\t'.join(['revision', 'user', 'date', 'systematic_id', 'primary_name', 'feature_type', 'qualifier_type', 'added_or_removed', 'value' ]) + '\n')
            out.write('\n'.join(sorted(qualifiers_output)))




