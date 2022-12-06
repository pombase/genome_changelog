from custom_biopython import SeqRecord, SeqFeature, SeqIO, CustomSeqFeature
import pandas
from copy import deepcopy
import warnings

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

    for normal_feature in genome.features:
        feature = CustomSeqFeature.from_parent(normal_feature, genome)
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
    - locations_added: list of SeqFeatures in new_genome_dict that have been added or their sequence has changed (sequence, not feature location since updates in genome sequence may shift multiple features, even the annotation is identical).
    - locations_removed: list of SeqFeatures in old_genome_dict that have been removed or their sequence has changed (sequence, not feature location since updates in genome sequence may shift multiple features, even the annotation is identical).
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

                ## First - changes to the feature sequence
                old_sequences = [f.feature_sequence for f in old_features]
                new_sequences = [f.feature_sequence for f in new_features]

                locations_added += [f for f in new_features if f.feature_sequence not in old_sequences]
                locations_removed += [f for f in old_features if f.feature_sequence not in new_sequences]

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

def make_synonym_dict(gene_ids_file):
    data = pandas.read_csv(gene_ids_file,sep='\t',na_filter=False)
    synonyms = dict()
    for i,row in data.iterrows():
        for synonym in [row['primary_name']] + row['synonyms'].split(','):
            if len(synonym):
                if synonym in synonyms:
                    synonyms[synonym].append(row['systematic_id'])
                else:
                    synonyms[synonym] = [row['systematic_id']]

    return synonyms


def read_pombe_genome(file_name, format, gene_ids_file, all_systematic_ids_ever, known_exceptions) -> SeqRecord:
    """
    Runs SeqIO.read on the file (with tweaks, see custom_biopython file), and then
    for features that do not have a systematic_id qualifier:
        * If they have a `\gene` qualifier with a current systematic id, add a systematic_id
        qualifier with that value.
        * If not, see if a `\gene` qualifier starting with `SP` is a UNIQUE synonym of a systematic_id, then use that.
        * If neither, see if there is a single `\gene` qualifier that starts by `SP`, this may mean that
        it corresponds to a systematic_id that was deleted.
    """
    known_exceptions_tsv = pandas.read_csv(known_exceptions,sep='\t')
    known_exception_dict = dict()
    for i,row in known_exceptions_tsv.iterrows():
        known_exception_dict[frozenset(row['gene_qualifiers'].split(','))] = row['change_to']

    with open(all_systematic_ids_ever) as f:
        valid_ids = set([line.rstrip('\n') for line in f])

    synonyms = make_synonym_dict(gene_ids_file)

    # Avoids some encoding errors
    with open(file_name, errors='replace') as ins:
        contig = SeqIO.read(ins,format)

    # Features to be added to the contig file (duplicated systematic_ids)
    extra_features = list()
    for feature in contig.features:
        if 'systematic_id' not in feature.qualifiers and 'gene' in feature.qualifiers:

            # The gene qualifier contains a current systematic_id (we use a set because sometimes there are duplicated qualifiers)
            systematic_ids = list(set(value for value in feature.qualifiers['gene'] if value in valid_ids))
            if len(systematic_ids) > 0:
                original_systematic_ids = systematic_ids
                # Substitute all values by their synonym + keep unique values only
                systematic_ids = list(set(synonyms[i][0] if i in synonyms else i for i in systematic_ids))

                for i in systematic_ids:
                    if i in synonyms and len(synonyms[i])> 1:
                        # Send a warning if we are using a synonym with multiple possibilities
                        warnings.warn(f'using a synonym with more than one possibility:\n{feature}')

                if len(systematic_ids) > 1:
                    if frozenset(original_systematic_ids) in known_exception_dict:
                        value = known_exception_dict[frozenset(original_systematic_ids)]
                        if value == 'duplicate':
                            original_systematic_ids = value[:1]
                            for systematic_id in value[1:]:
                                copied_feature = deepcopy(feature)
                                copied_feature.qualifiers['systematic_id'] = [systematic_id]
                                extra_features.append(copied_feature)
                        else:
                            systematic_ids = value
                    else:
                        raise ValueError('\\gene qualifier contains more than one systematic_id and not included in known_exceptions', systematic_ids)
                feature.qualifiers['systematic_id'] = systematic_ids
                continue

            # The gene qualifier contains a value starting with SP which is a unique synonym of a current systematic_id
            systematic_ids = set([value for value in feature.qualifiers['gene'] if (value.startswith('SP') and value in synonyms)])
            if len(systematic_ids) > 0:
                if len(systematic_ids) > 1:
                    print('skipped feature with gene qualifiers',feature.qualifiers['gene'],'because it\'s a synonym of two current systematic_ids', [synonyms[i] for i in feature.qualifiers['gene'] if i in synonyms])
                else:
                    feature.qualifiers['systematic_id'] = systematic_ids
                continue


            for value in feature.qualifiers['gene']:
                if value.startswith('SP'):
                    pass
                    # print('a value in \gene qualifier starts with SP but was never a systematic_id, skipping -> ', value)
                    # for qualifier_type in feature.qualifiers:
                    #     print('   ',qualifier_type,'-',feature.qualifiers[qualifier_type])

    contig.features.extend(extra_features)

    return contig
