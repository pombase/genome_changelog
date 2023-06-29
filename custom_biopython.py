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
    length_parts = re.match('(\d+)([^\d]+)',text).groups()
    assert len(length_parts) == 2, "Invalid sequence length string %r" % text
    assert length_parts[1].upper().strip() in ["BP", "BP.", "AA", "AA.", "P.", "B P."]
    consumer.size(length_parts[0])

def permissive_first_line_scanner(self, consumer, line):
    assert line[: self.HEADER_WIDTH].rstrip() == "ID"
    if line[self.HEADER_WIDTH :].count(";") == 6:
        # Looks like the semi colon separated style introduced in 2006
        self._feed_first_line_new(consumer, line)
    elif line[self.HEADER_WIDTH :].count(";") == 3:
        if line.rstrip().endswith(" SQ"):
            # EMBL-bank patent data
            self._feed_first_line_patents(consumer, line)
        else:
            # Looks like the pre 2006 style
            self._feed_first_line_old(consumer, line)
    elif line[self.HEADER_WIDTH :].count(";") == 2:
        # Looks like KIKO patent data
        self._feed_first_line_patents_kipo(consumer, line)
    # Special case for mating type region file
    elif 'XaviM' in line:
        self._feed_first_line(consumer, line + "; SV 1; linear; genomic DNA; STD; FUN; 19433 BP.")
    else:
        raise ValueError("Did not recognise the ID line layout:\n" + line)

EmblScanner._feed_seq_length = permissive_seq_length_scanner
EmblScanner._feed_first_line = permissive_first_line_scanner



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