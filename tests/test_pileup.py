import pysam
from mergedups.pileup import PileupElement, Pileup

__author__ = 'dankle'

my_read = pysam.AlignedSegment()
my_read.query_sequence="AGCT"
my_read.query_qualities = pysam.fromQualityString("FGHI")
my_read.reference_id = 0
my_read.reference_start = 32
my_read.mpos = 250
my_read.is_reverse = True
my_read.tags = (("NM", 1), ("RG", "L1"))

# my_read2 is truncated by one base
my_read2 = pysam.AlignedSegment()
my_read2.query_sequence="GGC"
my_read2.query_qualities = pysam.fromQualityString("FGH")

my_read3 = pysam.AlignedSegment()
my_read3.query_sequence="AACT"
my_read3.query_qualities = pysam.fromQualityString("FGHI")

my_read4 = pysam.AlignedSegment()
my_read4.query_sequence="GGCT"
my_read4.query_qualities = pysam.fromQualityString("FGHI")


def test_pileup_element():
    pe = PileupElement('A', 40)
    assert pe.phredqual() == 'I'
    assert pe.base == 'A'

    pe = PileupElement('G', None)
    assert pe.phredqual() == None
    assert pe.base == 'G'

    pe = PileupElement('T', phred='I')
    assert pe.qual == 40

def test_pileup():
    pe = PileupElement('A', 40)
    pe2 = PileupElement('A', 20)
    pe3 = PileupElement('G', None)
    pileup = Pileup([pe, pe2, pe3])
    assert pileup.bases() == ['A', 'A', 'G']
    assert pileup.quals() == [40, 20, None]
    assert pileup.most_common_base() == 'A'

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
    assert [len(p.bases()) for p in pileups] == [4,4,4,3]


