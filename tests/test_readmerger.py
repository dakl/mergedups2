import tempfile
import pysam
from mergedups.pileup import Pileup
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
my_read.is_paired = True
my_read.is_proper_pair = True

# my_read2 is truncated by one base
my_read2 = pysam.AlignedSegment()
my_read2.query_sequence="GGC"
my_read2.query_qualities = pysam.fromQualityString("FGH")
my_read2.is_paired = True
my_read2.is_proper_pair = True
my_read2.tags = (("NM", 1), ("RG", "L1"))

my_read3 = pysam.AlignedSegment()
my_read3.query_sequence="AACT"
my_read3.query_qualities = pysam.fromQualityString("FGHI")
my_read3.is_paired = True
my_read3.is_proper_pair = True
my_read3.tags = (("NM", 1), ("RG", "L1"))

my_read4 = pysam.AlignedSegment()
my_read4.query_sequence="GGCT"
my_read4.query_qualities = pysam.fromQualityString("FGHI")

header = {'HD': {'VN': '1.0'},
          'SQ': [{'LN': 1575, 'SN': 'chr1'},
                 {'LN': 1584, 'SN': 'chr2'}]}

a = pysam.AlignedSegment()
a.query_name = "read_28833_29006_6945"
a.query_sequence = "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
a.flag = 99
a.reference_id = 0
a.reference_start = 32
a.mapping_quality = 20
a.cigar = ((0, 10), (2, 1), (0, 25))
a.next_reference_id = 0
a.next_reference_start = 199
a.template_length = 167
a.query_qualities = pysam.fromQualityString("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF")
a.tags = (("NM", 1),
          ("RG", "L1"))
a.is_paired = True
a.is_proper_pair = True

b = pysam.AlignedSegment()
b.query_name = "read_28833_29006_6944"
b.query_sequence="GGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
b.flag = 99
b.reference_id = 0
b.reference_start = 32
b.mapping_quality = 20
b.cigar = ((0, 10), (2, 1), (0, 25))
b.next_reference_id = 0
b.next_reference_start = 199
b.template_length = 167
b.query_qualities = pysam.fromQualityString("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF")
b.tags = (("NM", 1),
          ("RG", "L1"))
b.is_paired = True
b.is_proper_pair = True

def test_key():
    expected_key = "0_32_250_True_L1"
    key = ReadMerger.key(my_read)
    assert key == expected_key

def test_merge():
    merged_pileup_elements = []
    """ :type: list[PileupElement] """
    pileups = []
    """ :type: list[Pileup] """
    reads = [my_read, my_read2, my_read3, my_read4]

    for k in range(0, 4):
        pile = Pileup.from_reads(reads, k)
        pileups.append(pile)
        merged_pileup_elements.append(pile.merge(max_qual=45, fraction_agree=0.75))

    # test that ties gives N
    assert merged_pileup_elements[0].base == 'N'  # tie, 2xA 2xG
    # test that majority vote works
    assert merged_pileup_elements[1].base == 'G'  # G wins from 3xG 1xA

    # test all 4 bases that were merged
    assert "".join([pe.base for pe in merged_pileup_elements]) == 'NGCT'
    # test that qualities were properly merged
    assert [pe.qual for pe in merged_pileup_elements] == [2, 45, 45, 45]
    # test that the size of each pileup works
    assert [len(p.bases()) for p in pileups] == [4, 4, 4, 3]

    assert ReadMerger.key(a) == ReadMerger.key(b)
    assert Pileup.from_reads([a, b], 0).merge(45, .75).base == 'N'
    assert Pileup.from_reads([a, b], 0).merge(45, .75).qual == 2
    assert Pileup.from_reads([a, b], 1).merge(45, .75).base == 'G'

    tmpfilename = tempfile.mkstemp('.bam')[1]
    outfile = pysam.AlignmentFile(tmpfilename, "wb", header=header)
    mf = "/tmp/bla.metrics.txt"
    readmerger = ReadMerger([a, b], outfile, mf, 0.75, 45, 1)
    exitcode = readmerger.do_work()
    assert exitcode == 0
    metrics = readmerger.metrics.dict()
    assert 'L1' in metrics
    # assert metrics['L1'] == {}
    #['L1']['READ_PAIR_DUPLICATES'] == 2

