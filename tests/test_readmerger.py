import pysam
from mergedups.readmerger import ReadMerger

__author__ = 'dankle'

my_read = pysam.AlignedSegment()
my_read.query_sequence="AGCT"
my_read.query_qualities = pysam.fromQualityString("FGHI")
my_read.reference_id = 0
my_read.reference_start = 32
my_read.mpos = 250
my_read.is_reverse = True
my_read.tags = (("NM", 1), ("RG", "L1"))

def test_key():
    expected_key = "0_32_250_True_L1"
    key = ReadMerger.key(my_read)
    assert key == expected_key
