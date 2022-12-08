"""
Some overwritten biopython functions
"""

from Bio.GenBank.Scanner import EmblScanner
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
# Included here to make sure EmblScanner._feed_seq_length is overwritten
from Bio import SeqIO
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


def features_are_equal(self: SeqFeature, other):
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

class CustomSeqFeature(SeqFeature):

    def __init__(self, location=None, type="", location_operator="", strand=None, id="<unknown id>", qualifiers=None, sub_features=None, ref=None, ref_db=None, reference_sequence: SeqRecord =None):
        super().__init__(location, type, location_operator, strand, id, qualifiers, sub_features, ref, ref_db)
        self.reference_sequence = reference_sequence
        self.feature_sequence = self.extract(self.reference_sequence).seq

    def __eq__(self, other):
        return super().__eq__(other) and (self.reference_sequence.seq == other.reference_sequence.seq)

    @classmethod
    def from_parent(cls, parent_instance: SeqFeature, reference_sequence: SeqRecord):
        return cls(parent_instance.location, parent_instance.type, parent_instance.location_operator, parent_instance.strand, parent_instance.id, parent_instance.qualifiers, None, parent_instance.ref, parent_instance.ref_db, reference_sequence)