"""
Containers for merge duplicate metrics, such as nubmer of reads,
number of dups, pct dups etc
"""

__author__ = 'dankle'

class ReadMergerMetrics(object):
    """
    contains metrics from a bam file. Holds a LibraryMetrics instance
    for each RG tag found in the bam file
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

    def update_pct_duplicated(self):
        """
        updated the estimated percent duplication
        :return:
        """
        for readgroup in self.libraries:
            self.libraries[readgroup].update_pct_duplicated()

    def dict(self):
        """
        Return a dictionary representation of the ReadMergerMetrics object
        """
        as_dict = {}
        for lib in self.libraries:
            as_dict[lib] = self.libraries[lib].dict()
        return as_dict

class LibraryMetrics(object):
    """
    Container for metrics from a single library
    """
    def __init__(self):
        self.stats = {}
        self.unpaired_reads_examined = 0
        self.read_pairs_examined = 0
        self.unmapped_reads = 0
        self.unpaired_read_duplicates = 0
        self.read_pair_duplicates = 0
        self.read_pair_optical_duplicates = 0
        self.percent_duplication = 0.0
        self.estimated_library_size = 0

    def update_metrics(self, read):
        """
        update metrics based on the read
        :type read: pysam.AlignedSegment
        :return: None
        """
        if read.is_unmapped:
            self.unmapped_reads += 1
        elif (not read.is_paired) or \
                read.mate_is_unmapped:  # read is single end or mate is unmapped
            self.unpaired_reads_examined += 1
            if read.is_duplicate:
                self.unpaired_read_duplicates += 1
        else:  # read is PE and both mates are mapped
            if read.is_read1:  # only update metrics for reads that are the first in a pair.
                self.read_pairs_examined += 1
                if read.is_duplicate:
                    self.read_pair_duplicates += 1

    def update_pct_duplicated(self):
        """
        updated the estimated percent duplication
        :return:
        """
        if self.read_pairs_examined+self.unpaired_reads_examined > 0:
            numer = 0.0 + self.read_pair_duplicates + self.unpaired_read_duplicates
            denom = 0.0 + self.read_pairs_examined + self.unpaired_reads_examined
            self.percent_duplication = numer / denom

    def dict(self):
        """
        Return a dictionary representation of the LibraryMetrics object
        :return:
        """
        return {'UNPAIRED_READS_EXAMINED':self.unpaired_reads_examined,
                'READ_PAIRS_EXAMINED': self.read_pairs_examined,
                'UNMAPPED_READS': self.unmapped_reads,
                'UNPAIRED_READ_DUPLICATES': self.unpaired_read_duplicates,
                'READ_PAIR_DUPLICATES': self.read_pair_duplicates,
                'READ_PAIR_OPTICAL_DUPLICATES': self.read_pair_optical_duplicates,
                'PERCENT_DUPLICATION': self.percent_duplication,
                'ESTIMATED_LIBRARY_SIZE': self.estimated_library_size}
