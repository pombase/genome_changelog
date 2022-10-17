#%%
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import glob
from Bio.SeqIO.InsdcIO import _insdc_location_string as format_location

def features_are_equal(self, other):
    if type(self) != type(other):
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

## Override print location
FeatureLocation.__str__ = lambda x: format_location(x, None)

def build_seqfeature_dict(genome: SeqRecord):

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

#TODO proper sorting
all_files = sorted(glob.glob('revision_files/chromosome*.txt'))
new_genome_dict = build_seqfeature_dict(SeqIO.read(all_files.pop(),'embl'))
old_genome_dict = build_seqfeature_dict(SeqIO.read(all_files.pop(),'embl'))

with open('revisions_chromosome1.txt') as f:
    revisions = f.read().splitlines()[::-1]


print('revision','user', 'date','systematic_id', 'feature_type', 'change', 'old_value', 'new_value',sep='\t')

def format_location_change(revision, features, label):
    out_str = ''
    for f in features:
        print(*revision, systematic_id,*format_location_change(old_genome_dict[systematic_id], 'removed'))
    return out_str




while new_genome_dict:
    revision = revisions.pop().split()
    locations_added = list()
    locations_removed = list()

    # Only in the new
    for systematic_id in set(new_genome_dict.keys()) - set(old_genome_dict.keys()):
        locations_added += list(new_genome_dict[systematic_id].values())

    # Only in the old
    for systematic_id in set(old_genome_dict.keys()) - set(new_genome_dict.keys()):
        locations_removed += list(old_genome_dict[systematic_id].values())

    # The shared ones
    for systematic_id in set(old_genome_dict.keys()).intersection(set(new_genome_dict.keys())):

        if old_genome_dict[systematic_id] != new_genome_dict[systematic_id]:
            new_annotation = new_genome_dict[systematic_id]
            old_annotation = old_genome_dict[systematic_id]

            # Only in the new
            for feature_type in set(new_annotation.keys()) - set(old_annotation.keys()):
                locations_added += new_genome_dict[systematic_id][feature_type]

            # Only in the old
            for feature_type in set(old_annotation.keys()) - set(new_annotation.keys()):
                locations_removed += old_genome_dict[systematic_id][feature_type]

            # The shared ones
            for feature_type in set(old_annotation.keys()).intersection(new_annotation.keys()):
                old_locations = [f.location for f in old_annotation[feature_type]]
                new_locations = [f.location for f in new_annotation[feature_type]]

                locations_added += [f for f in sum(list(new_genome_dict[systematic_id].values()),[]) if f.location not in old_locations]
                locations_removed += [f for f in sum(list(old_genome_dict[systematic_id].values()),[]) if f.location not in new_locations]

                old_qualifiers = set()
                new_qualifiers = set()

                for f in old_annotation[feature_type]:
                    old_qualifiers.update( set([(key, tuple(value)) for key, value in f.qualifiers.items()]))
                for f in new_annotation[feature_type]:
                    new_qualifiers.update( set([(key, tuple(value)) for key, value in f.qualifiers.items()]))

                for key, value in new_qualifiers - old_qualifiers:
                    print(*revision, systematic_id, feature_type, 'qualifier:added',  f'{key}:{value}' , sep='\t')
                for key, value in old_qualifiers - new_qualifiers:
                    print(*revision, systematic_id, feature_type, 'qualifier:added',  f'{key}:{value}' , sep='\t')

    old_genome_dict = new_genome_dict
    if len(all_files):
        new_genome_dict = build_seqfeature_dict(SeqIO.read(all_files.pop(),'embl'))
    else:
        break



