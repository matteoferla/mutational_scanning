#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Written for python 3, not tested under 2.
"""
QQC in Python. Under construction. See the Matlab files!
"""
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "25/10/2016"

N = "\n"
T = "\t"
# N = "<br/>
import re
from warnings import warn


class Trace:
    bases = ('A', 'T', 'G', 'C')

    def __init__(self, record):
        for ni, N in enumerate(record.annotations['abif_raw']['FWO_1'].upper()):
            setattr(self, N, record.annotations['abif_raw']['DATA{0}'.format(9 + ni)])
        self.peak_index = record.annotations['abif_raw']['PLOC1']
        self.peak_id = record.annotations['abif_raw']['PBAS1']

    @classmethod
    def from_filename(cls, file):
        """
        Class method to return a Trace instance, but from a filename.
        >>> data=Trace.from_filename('sample.ab1')
        :param file: a file name
        :return: Trace instance
        """
        handle = open(file, 'rb')
        record = SeqIO.read(handle, "abi")
        return cls(record)

    def find_peak(self, target_seq, strict=True):
        """
        Find the index of the first base of a given sequence in the trace.
        :param target_seq: a string of bases.
        :param strict: boolean, dies if multiple hits.
        :return: the index of peak_id or peak_index of the first base.
        """
        ## sanitise
        target_seq = re.sub('[^ATGC]', '', target_seq.upper())
        targetdex = re.findall(target_seq, self.peak_id)
        if not targetdex:
            raise ValueError('target_seq not found!')
        elif len(targetdex) > 1:
            if strict:
                raise ValueError('Ambiguous sequence given!')
            else:
                warn('Ambiguous sequence given!')
        return self.peak_id.find(target_seq)

    def find_peak_after(self, target_seq, strict=True):
        """
        Same as find_peak, but give the index of the base after the given sequence.
        :param target_seq: a string of bases.
        :param strict: boolean, dies if multiple hits.
        :return: the index of peak_id or peak_index of the next (last+1) base along
        """
        target_seq = re.sub('[^ATGC]', '', target_seq.upper())
        return self.find_peak(target_seq, strict)+len(target_seq)+1

    def get_intensities(self, index,wobble=0.20):
        """
        Get intensities at given peak index.
        :param index: a peak index
        :param wobble: window to look in (default .2)
        :return: dictionary
        """
        span = round(len(self.A)/len(self.peak_index))
        doubleindex=self.peak_index[index]
        return {base: max([getattr(self,base)[doubleindex+i] for i in range(-round(span*wobble), round(span*wobble))]) for base in self.bases}

    def QQC(self,location, *args, **kwargs):
        if isinstance(location, str):
            location = self.find_peak_after(location)
            ## go..
        QQC([self.get_intensities(location),self.get_intensities(location+1),self.get_intensities(location+2)],*args, **kwargs)


### Considering...
# I thought this would be the best way, but it is overkill.
class CodonProbs:
    bases = ('A', 'T', 'G', 'C')
    def __init__(self,codon):
        self._data=codon

    @classmethod
    def from_zeros(cls):
        return cls([{b:0 for b in cls.bases} for i in range(3)])

    def __iter__(self):
        for i in range(3):
            for b in self.bases:
                yield self._data[i][b]

    def __add__(self, other):
        return CodonProbs([{b: self._data[i][b] + other for b in self.bases} for i in range(3)])

    def __abs__(self):
        return CodonProbs([{b: abs(self._data[i][b]) for b in self.bases} for i in range(3)])


    def clone(self):
        return CodonProbs(self._data)

    def bsxfun(self, fun):
        new=CodonProbs.from_zeros()
        for i in range(3):
            for b in self.bases:
                new._data[i][b]=fun(self._data[i][b])
        return new
#### end...

class QQC:
    def __init__(self,peak_int,scheme='NNK'):
        self.scheme = scheme
        scheme_pred=[{b:.25 for b in 'ATGC'},{b:.25 for b in 'ATGC'},{'A':0,'T':.50, 'G':.50,'C':0}]
        # make peakfreqs a fraction of one
        codon_peak_freq=[{b:p[b]/sum(p.values()) for b in 'ATGC'} for p in peak_int]
        # calculate deviation per position
        deviation=[sum([scheme_pred[i][b] - abs(scheme_pred[i][b] - codon_peak_freq[i][b]) for b in 'ATGC']) for i in range(3)]
        weight_total=sum([scheme_pred[i][b]>0 for b in 'ATGC' for i in range(3)])
        weights = [sum([scheme_pred[i][b]>0 for b in 'ATGC'])/weight_total for i in range(3)]
        wsum=sum([deviation[i] * weights[i] for i in range(3)])
        undiverse=[{'A':1,'T':0, 'G':0,'C':0},{'A':1,'T':0, 'G':0,'C':0},{'A':1,'T':0, 'G':0,'C':0}]
        worse = [sum([scheme_pred[i][b] - abs(scheme_pred[i][b] - undiverse[i][b]) for b in 'ATGC']) for i in range(3)]
        wmin=sum([worse[i] * weights[i] for i in range(3)])
        self.Qpool=(wsum+abs(wmin))/(1+abs(wmin))

    @staticmethod
    def from_trace(trace, location, *args, **kwargs):
        """
        An alternate route. QQC instance comes from Trace(...).QQC()
        :param trace: Trace instance
        :param location: location of first peak of three
        :return: QQC instance
        """
        if not isinstance(trace, Trace):
            trace = Trace(trace)
        return trace.QQC(location, *args, **kwargs)










if __name__ == "__main__":
    # the documentation is shocking
    # http://biopython.org/DIST/docs/api/Bio.SeqIO.AbiIO-module.html
    from Bio import SeqIO
    import json

    file = "example data/ACE-AA-088-01-55Â°C-BM3-A82_19C-T7-T7minus1.ab1"
    x = Trace.from_filename(file)
    print(x.QQC('CGT GAT TTT'))
    # To Do...
    # get base frequency at given peak
    # QQC class
    # given input CGT GAT TTT NNK
    # calculate ratios
    #
