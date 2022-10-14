#%%
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
import glob
def features_are_equal(self, other):

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

    out_dict = dict()

    for feature in genome.features:
        feature: SeqFeature
        if 'systematic_id' not in feature.qualifiers:
            continue
        gene_id = feature.qualifiers['systematic_id'][0]
        if gene_id not in out_dict:
            out_dict[gene_id] = dict()
        feature_type = feature.type
        if feature_type in ['intron', 'misc_feature']:
            continue
        if feature_type in out_dict[gene_id]:
            raise ValueError(
                f'several features of {feature_type} for {gene_id}')

        out_dict[gene_id][feature_type] = feature

    return out_dict



#TODO proper sorting
all_files = sorted(glob.glob('revision_files/chromosome*.txt'))
new_features = build_seqfeature_dict(SeqIO.read(all_files.pop(),'embl'))
old_features = build_seqfeature_dict(SeqIO.read(all_files.pop(),'embl'))
revisions = [f.replace('chromosome1_','').replace('.txt','') for f in all_files]

print('revision','systematic_id', 'change' , 'feature_type', 'old_feature_location', 'new_feature_location',sep='\t')

while new_features:
    revision = revisions.pop()
    for systematic_id in set(list(new_features.keys()) + list(old_features.keys())):
        if systematic_id not in new_features:
            pass
            print(revision, systematic_id, 'removed', sep='\t')
        elif systematic_id not in old_features:
            pass
            print(revision, systematic_id, 'added', sep='\t')
        else:
            if old_features[systematic_id] != new_features[systematic_id]:
                new_annotation = new_features[systematic_id]
                old_annotation = old_features[systematic_id]
                for feature_type in set(list(new_annotation.keys()) + list(old_annotation.keys())):
                    if feature_type not in new_annotation:
                        pass
                        print(revision, systematic_id, 'removed', feature_type, sep='\t')
                    elif feature_type not in old_annotation:
                        pass
                        print(revision, systematic_id, 'added', feature_type, sep='\t')
                    else:
                        # print(revision, k, f'changed from {old_annotation[k]} to {new_annotation[k]}')
                        old_feature = old_annotation[feature_type]
                        new_feature = new_annotation[feature_type]
                        if old_feature.qualifiers != new_feature.qualifiers:
                            pass
                            print(revision, systematic_id, 'qualifiers changed', feature_type, sep='\t')
                        if old_feature.location != new_feature.location:
                            print(revision, systematic_id, 'modified' , feature_type, old_feature.location, new_feature.location, sep='\t')

    old_features = new_features
    if len(all_files):
        try:
            new_features = build_seqfeature_dict(SeqIO.read(all_files.pop(),'embl'))
        except ValueError as e:
            print(revision, e, sep='\t')
    else:
        break



