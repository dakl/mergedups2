import pysam
import logging
from datetime import datetime
from mergedups.pileup import Pileup
from mergedups.metrics import ReadMergerMetrics

__author__ = 'dankle'


class ReadMerger(object):
    """
    Class that merges `pysam.AlignedSegment`s.
    """

    def __init__(self, input_bam, output_bam,
                 fraction_agree, max_qual, reads_between_logs):
        """
        :type input_bam: pysam.AlignmentFile
        :type output_bam: pysam.AlignmentFile
        :type fraction_agree: float
        :type max_qual: int
        :type reads_between_logs: int
        :return: ReadMerger
        """
        self.input_bam = input_bam
        self.output_bam = output_bam
        self.fraction_agree = fraction_agree
        self.max_qual = max_qual
        self.reads_between_logs = reads_between_logs
        self.start_time = None
        self.last_time = None
        self.n_reads_traversed = 0
        self.last_n_reads_traversed = 0
        self.metrics = ReadMergerMetrics()
        self.last_chr = None
        self.last_pos = None
        self.name_translation_table = {}
        """ :type: dict[str,str] """

    def do_work(self):
        """
        Traverse the input_bam object and write merged reads to output_bam.
        :return: int
        """
        self.start_time = datetime.now()
        self.last_time = self.start_time
        reads = {}
        """:type: dict[str,list[pysam.AlignedSegment]]"""

        for read in self.input_bam:
            if not read.is_paired and not read.is_proper_pair:
                logging.debug("Skipping read {}".format(ReadMerger.key(read)))
                continue

            if self.last_chr is None or self.last_pos is None:
                self.last_chr = read.reference_id
                self.last_pos = read.pos

            key = ReadMerger.key(read)
            if key in reads:
                reads[key].append(read)
            else:
                reads[key] = [read]
                logging.debug("Added 1 read to reads {}".format(reads[key]))

            if read.pos > self.last_pos or read.reference_id != self.last_chr:
                logging.debug("Current chr:pos is {}:{}".format(read.reference_id, read.pos))
                self.work_on_reads(reads)
                self.last_pos = read.pos
                self.last_chr = read.reference_id

            self.n_reads_traversed += 1

            if self.n_reads_traversed % self.reads_between_logs == 0:
                self.log_traversal_status()

        self.work_on_reads(reads)

        self.log_final_status()
        self.metrics.update_pct_duplicated()
        logging.info(self.metrics.dict())
        return 0

    def work_on_reads(self, reads):
        """
        Check if there are any duplicates to be merged in reads
        :param reads: dict[str, list[pysam.AlignedSegment]]
        :return: None
        """
        keys = reads.keys()
        for key in keys:
            reads_to_merge = reads[key]
            for r in reads_to_merge:
                self.name_translation_table[r.qname] = reads_to_merge[0].qname
            logging.debug("Key is {}".format(key))
            logging.debug("Reads to merge are {}".format(reads_to_merge))
            merged_read = ReadMerger.merge(reads[key], self.max_qual, self.fraction_agree)

            # rename read so that pairs match in the final file
            if merged_read.qname in self.name_translation_table:
                merged_read.qname = self.name_translation_table[merged_read.qname]

            self.output_bam.write(merged_read)
            self.update_metrics(reads[key])
            del reads[key]



    def update_metrics(self, reads):
        """
        :type reads: list[pysam.AlignedSegment]
        :return: None
        """
        logging.debug("Updating metrics for {} reads".format(len(reads)))
        # reads 1..len(reads) are duplicates
        dup_reads = [reads[k] for k in range(1, len(reads))]
        for read in dup_reads:
            read.is_duplicate = True
            self.metrics.update_metrics(read)

        # of only a single read was found, it's not a duplicate
        self.metrics.update_metrics(reads[0])

    def log_traversal_status(self):
        """
        Log traversal status while running
        :return:
        """
        curr_time = datetime.now()
        time_since_start = curr_time - self.start_time
        time_since_start = time_since_start.seconds + time_since_start.microseconds/1E6
        timediff = curr_time - self.last_time
        timediff = timediff.seconds + timediff.microseconds/1E6
        if timediff == 0:  # don't divide by 0
            timediff = 1.0
        n_reads_since_last = self.n_reads_traversed - self.last_n_reads_traversed
        logging.info("Traversed {n} reads in {t} seconds at {x} reads per second".format(
            n=self.n_reads_traversed, t=time_since_start, x=n_reads_since_last/timediff
            ))
        self.last_time = curr_time
        self.last_n_reads_traversed = self.n_reads_traversed

    def log_final_status(self):
        """
        Log final traversal stats
        :return:
        """
        curr_time = datetime.now()
        time_since_start = curr_time - self.start_time
        time_since_start = time_since_start.seconds + time_since_start.microseconds/1E6
        logging.info(
            "Traversed total {n} reads in {t} seconds at an average rate of {x} reads per second".format(
                n=self.n_reads_traversed, t=time_since_start, x=self.n_reads_traversed/time_since_start
            ))
        self.last_time = curr_time
        self.last_n_reads_traversed = self.n_reads_traversed

    @staticmethod
    def merge(reads, max_qual=46, fraction_agree=0.75):
        """
        :type reads: list[pysam.AlignedSegment]
        :type max_qual: int
        :type fraction_agree: float
        :return: pysam.AlignedSegment
        """
        logging.debug("Merging {} read(s).".format(len(reads)))
        if len(reads) == 1:
            return reads[0]
        merged = reads[0]
        max_read_len = max([len(r.seq) for r in reads])
        new_seq = []
        new_quals = []
        for k in range(max_read_len):
            pileup = Pileup.from_reads(reads, k)
            merged_pileup_element = pileup.merge(max_qual, fraction_agree)
            new_seq.append(merged_pileup_element.base)
            new_quals.append(merged_pileup_element.phredqual())

        logging.debug("new seq is {}".format("".join(new_seq)))
        logging.debug("new qual is {}".format("".join(new_quals)))
        merged.seq = "".join(new_seq)
        merged.qual = "".join(new_quals)
        return merged

    @staticmethod
    def key(read):
        """
        Get a unique identifier read based on chr/pos/strand of the read and it's mate.
        If two reads have identical keys, they are PCR duplicates.
        :type read: pysam.AlignedSegment
        :return: str
        """

        # query_alignment_start is the index of the first
        # non-soft-clipped base in the sequence (0 for non-soft-clipped reads)
        readgroup = None
        if read.has_tag('RG'):
            readgroup = read.get_tag('RG')
        return "{chr}_{pos}_{matepos}_{strand}_{readgroup}".format(
            chr=read.reference_id,
            pos=read.pos-read.query_alignment_start,
            matepos=read.mpos,
            strand=read.is_reverse,
            readgroup=readgroup
        )

