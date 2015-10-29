import logging

__author__ = 'dankle'

class PileupElement(object):
    """
    Container for a single base with an associated numeric quality metric. The underlying read is
    optionally available as e.read and the associated offset is available in e.offset
    """
    def __init__(self, base, qual=None, phred=None, read=None, offset=None):
        self.base = base
        self.qual = None
        if qual:
            self.qual = qual
        if phred:
            self.qual = ord(phred) - 33
        self.read = read
        self.offset = offset

    def phredqual(self):
        if self.qual is not None:
            return chr(self.qual + 33)
        else:
            return None

class Pileup(object):
    def __init__(self, pileup_elements):
        self.elements = pileup_elements

    def bases(self):
        """
        Returns a list of the bases in the pileup elements
        :return: list[str]
        """
        return [e.base for e in self.elements]

    def quals(self):
        """
        Returns a list of base qualities for the pileup elements
        :return: list[int]
        """
        return [e.qual for e in self.elements]

    def most_common_base(self):
        """
        Returns the most frequently occurring base in the pileup as a string.
        If there are ties, the first of the two most frequently occurring elements is returned.
        :return: str
        """
        return max(set(self.bases()), key=self.bases().count)

    def merge(self, max_qual, fraction_agree):
        """
        :type max_qual: int
        :type fraction_agree: float
        :return: PileupElement
        """
        (base,qual_sum) = ('N',2)
        bases_set = set(self.bases())

        if len(bases_set) == 1:
            base = bases_set.pop()
            qual_sum = sum(self.quals())
            if qual_sum > max_qual:
                qual_sum = max_qual
        else:
            most_common_base = self.most_common_base()
            # which(bases == most_common_base) in python
            idx = [i for i,x in enumerate(self.bases()) if x == most_common_base]

            # if at least the given fraction of the reads support the most common base
            if len(idx) >= fraction_agree*len(self.bases()):
                qual_sum = sum([self.quals()[i] for i in idx])
                if qual_sum > max_qual:
                    qual_sum = max_qual
                base = most_common_base

            else:
                base = 'N'
                qual_sum = 2

        pileup_element = PileupElement(base, qual_sum)
        return pileup_element

    @staticmethod
    def from_reads(reads, pos):
        return Pileup([PileupElement(read.seq[pos], phred=read.qual[pos]) for read in reads if len(read.seq) > pos])
