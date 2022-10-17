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
        if feature_type in out_dict[gene_id]:
            if type(out_dict[gene_id][feature_type]) != list:
                out_dict[gene_id][feature_type] = [out_dict[gene_id][feature_type]]
            out_dict[gene_id][feature_type].append(feature)
        else:
            out_dict[gene_id][feature_type] = feature

    return out_dict

#TODO proper sorting
all_files = sorted(glob.glob('revision_files/chromosome*.txt'))
new_features = build_seqfeature_dict(SeqIO.read(all_files.pop(),'embl'))
old_features = build_seqfeature_dict(SeqIO.read(all_files.pop(),'embl'))

with open('revisions_chromosome1.txt') as f:
    revisions = f.read().splitlines()[::-1]


print('revision','user', 'date','systematic_id', 'feature_type', 'change', 'old_value', 'new_value',sep='\t')

def format_location_change(feature, label):
    if type(feature) == list:
        for f in feature:
            return f'list:{label}', f.location
    else:
        return label, feature.location



while new_features:
    revision = revisions.pop().split()
    for systematic_id in set(list(new_features.keys()) + list(old_features.keys())):
        if systematic_id not in new_features:
            for feature_type in old_features[systematic_id]:
                print(*revision, systematic_id,*format_location_change(old_features[systematic_id][feature_type], 'removed'))
        elif systematic_id not in old_features:
            for feature_type in new_features[systematic_id]:
                print(*revision, systematic_id,*format_location_change(new_features[systematic_id][feature_type], 'added'))
        else:
            if old_features[systematic_id] != new_features[systematic_id]:
                new_annotation = new_features[systematic_id]
                old_annotation = old_features[systematic_id]
                for feature_type in set(list(new_annotation.keys()) + list(old_annotation.keys())):
                    if feature_type not in new_annotation:
                        print(*revision, systematic_id,*format_location_change(old_annotation[feature_type], 'removed'))
                    elif feature_type not in old_annotation:
                        print(*revision, systematic_id,*format_location_change(new_annotation[feature_type], 'added'))
                    else:
                        old_feature = old_annotation[feature_type]
                        new_feature = new_annotation[feature_type]
                        print(*revision, systematic_id,*format_location_change(old_annotation[feature_type], 'removed'))
                        print(*revision, systematic_id,*format_location_change(new_annotation[feature_type], 'added'))
                        # old_qualifiers = set([old_feature.qualifiers]) if type(old_feature.qualifiers) != list else set([f.qualifiers])
                        if old_feature.qualifiers != new_feature.qualifiers:
                            for new_qualifier_name in new_feature.qualifiers:
                                if (new_qualifier_name not in old_feature.qualifiers) or (old_feature.qualifiers[new_qualifier_name] != new_feature.qualifiers[new_qualifier_name]):
                                    print(*revision, systematic_id, feature_type, 'qualifier:added',  f'{new_qualifier_name}:{new_feature.qualifiers[new_qualifier_name]}' , sep='\t')
                            for old_qualifier_name in old_feature.qualifiers:
                                if (old_qualifier_name not in new_feature.qualifiers) or (old_feature.qualifiers[old_qualifier_name] != new_feature.qualifiers[old_qualifier_name]):
                                    print(*revision, systematic_id, feature_type, 'qualifier:removed', f'{old_qualifier_name}:{old_feature.qualifiers[old_qualifier_name]}' , sep='\t')

    old_features = new_features
    if len(all_files):
        new_features = build_seqfeature_dict(SeqIO.read(all_files.pop(),'embl'))
    else:
        break



