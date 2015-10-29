"""
Containers for merge duplicate metrics, such as nubmer of reads, number of dups, pct dups etc
"""
import pysam

__author__ = 'dankle'

class ReadMergerMetrics(object):
    """
    contains metrics from a bam file
    """
    def __init__(self):
        self.libraries = {}
        """ :type : dict[str,LibraryMetrics] """

    def update_metrics(self, read):
        """
        update the metrics base on the given read
        :type read: pysam.AlignedSegment
        :return: Non
        """
        readgroup = read.get_tag('RG')
        if not readgroup in self.libraries:
            self.libraries[readgroup] = LibraryMetrics()
        self.libraries[readgroup].update_metrics(read)

    def dict(self):
        r = {}
        for lib in self.libraries:
            r[lib] = self.libraries[lib].dict()
        return r

class LibraryMetrics(object):
    """
    Container for metrics from a single library
    """
    def __init__(self):
        self.stats = {}
        self.UNPAIRED_READS_EXAMINED = 0
        self.READ_PAIRS_EXAMINED = 0
        self.UNMAPPED_READS = 0
        self.UNPAIRED_READ_DUPLICATES = 0
        self.READ_PAIR_DUPLICATES = 0
        self.READ_PAIR_OPTICAL_DUPLICATES = 0
        self.PERCENT_DUPLICATION = 0.0
        self.ESTIMATED_LIBRARY_SIZE = 0

    def update_metrics(self, read):
        """
        update metrics based on the read
        :type read: pysam.AlignedSegment
        :return: None
        """
        if read.is_unmapped:
            self.UNMAPPED_READS += 1
        elif (not read.is_paired) or read.mate_is_unmapped:  ## read is single end or mate is unmapped
            self.UNPAIRED_READS_EXAMINED
            if read.is_duplicate:
                self.UNPAIRED_READ_DUPLICATES += 1
        else: ## read is PE and both mates are mapped
            if read.is_read1: ## only update metrics for reads that are the first in a pair.
                self.READ_PAIRS_EXAMINED += 1
                if read.is_duplicate:
                    self.READ_PAIR_DUPLICATES += 1

        if self.READ_PAIRS_EXAMINED+self.UNPAIRED_READS_EXAMINED > 0:
            self.PERCENT_DUPLICATION = (0.0 + self.READ_PAIR_DUPLICATES + self.UNPAIRED_READ_DUPLICATES) / (0.0 + self.READ_PAIRS_EXAMINED+self.UNPAIRED_READS_EXAMINED)

    def dict(self):
        return {'UNPAIRED_READS_EXAMINED':self.UNPAIRED_READS_EXAMINED,
                'READ_PAIRS_EXAMINED': self.READ_PAIRS_EXAMINED,
                'UNMAPPED_READS': self.UNMAPPED_READS,
                'UNPAIRED_READ_DUPLICATES': self.UNPAIRED_READ_DUPLICATES,
                'READ_PAIR_DUPLICATES': self.READ_PAIR_DUPLICATES,
                'READ_PAIR_OPTICAL_DUPLICATES': self.READ_PAIR_OPTICAL_DUPLICATES,
                'PERCENT_DUPLICATION': self.PERCENT_DUPLICATION,
                'ESTIMATED_LIBRARY_SIZE': self.ESTIMATED_LIBRARY_SIZE
        }