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
        return self.find_peak(target_seq, strict) + len(target_seq) + 1

    def get_intensities(self, index, wobble=0.20):
        """
        Get intensities at given peak index.
        :param index: a peak index
        :param wobble: window to look in (default .2)
        :return: dictionary
        """
        span = round(len(self.A) / len(self.peak_index))
        doubleindex = self.peak_index[index]
        return {
        base: max([getattr(self, base)[doubleindex + i] for i in range(-round(span * wobble), round(span * wobble))])
        for base in self.bases}

    def QQC(self, location, *args, **kwargs):
        if isinstance(location, str):
            location = self.find_peak_after(location)
            ## go..
        return QQC(
            [self.get_intensities(location), self.get_intensities(location + 1), self.get_intensities(location + 2)],
            *args, **kwargs)


### Considering...
# I thought this would be the best way, but it is overkill.
class CodonProbs:
    bases = ('A', 'T', 'G', 'C')

    def __init__(self, codon):
        self._data = codon

    @classmethod
    def from_zeros(cls):
        return cls([{b: 0 for b in cls.bases} for i in range(3)])

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
        new = CodonProbs.from_zeros()
        for i in range(3):
            for b in self.bases:
                new._data[i][b] = fun(self._data[i][b])
        return new


#### end...

class QQC:
    def __init__(self, peak_int, scheme='NNK'):
        self.scheme = scheme
        self.scheme_mix = QQC.scheme_maker(scheme)
        self.scheme_pred = [
            {b: sum([self.scheme_mix[p][0] * self.scheme_mix[p][1][i][b] for p in range(len(self.scheme_mix))]) for b in
             'ATGC'} for i in range(3)]
        # make peakfreqs a fraction of one
        codon_peak_freq = [{b: p[b] / sum(p.values()) for b in 'ATGC'} for p in peak_int]
        # calculate deviation per position
        deviation = [sum([self.scheme_pred[i][b] - abs(self.scheme_pred[i][b] - codon_peak_freq[i][b]) for b in 'ATGC'])
                     for i in range(3)]
        # calculate the weights
        weight_total = sum([self.scheme_pred[i][b] > 0 for b in 'ATGC' for i in range(3)])
        weights = [sum([self.scheme_pred[i][b] > 0 for b in 'ATGC']) / weight_total for i in range(3)]
        # weighted sum of deviations
        wsum = sum([deviation[i] * weights[i] for i in range(3)])
        # calculated the worse case scenario, an undiverse codon
        undiverse = [{'A': 1, 'T': 0, 'G': 0, 'C': 0}, {'A': 1, 'T': 0, 'G': 0, 'C': 0},
                     {'A': 1, 'T': 0, 'G': 0, 'C': 0}]
        worse = [sum([self.scheme_pred[i][b] - abs(self.scheme_pred[i][b] - undiverse[i][b]) for b in 'ATGC']) for i in
                 range(3)]
        wmin = sum([worse[i] * weights[i] for i in range(3)])
        # calculate Qpool
        self.Qpool = (wsum + abs(wmin)) / (1 + abs(wmin))

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

    @staticmethod
    def scheme_maker(scheme):
        """
        converts a name of a scheme to a probability set.
        The input is a string (e.g. NNK). Except some are complicated.
        So it accepts one or more codons separated by a space and optionally prefixed with a interger number.
        >>> QQC.scheme_maker('12NDT 6VHA 1TGG 1ATG')
        :param scheme: str of codon proportion.
        :return: a probability list-dictionary thinggy
        """
        # check first if reserved word.
        scheme=scheme.upper().replace('-trick'.upper(),'')
        if scheme == 'Tang'.upper():
            scheme='12NDT 6VHA 1TGG 1ATG'
        elif scheme.lower() == '19c':
            scheme =''
        elif scheme.lower() == '20c':
            scheme =''
        elif scheme.lower() == '21c':
            scheme =''
        elif scheme.lower() == '22c':
            scheme ='1NDT 9VHG 1TGG'
        else:
            pass #there seems no need to store a boolean?

        degeneracy = {'N': 'ATGC',
                      'A': 'A',
                      'T': 'T',
                      'G': 'G',
                      'C': 'C',
                      'W': 'AT',
                      'S': 'GC',
                      'K': 'GT',
                      'M': 'AC',
                      'R': 'AG',
                      'Y': 'TC',
                      'B': 'TGC',
                      'V': 'AGC',
                      'H': 'ATC',
                      'D': 'ATG'}
        schemelist = scheme.split()
        proportions = []
        codons = []
        for m in schemelist:
            (prop, codon) = re.match('(\d{0,2})(\w{3})', m).groups()
            if prop:
                proportions.append(int(prop))
            else:
                proportions.append(1)
            codons.append(codon)
        freq = [p / sum(proportions) for p in proportions]
        codonmix = []
        for codon in codons:
            pred = []
            for i in range(3):
                basedex = {}
                for b in 'ATGC':
                    if b in degeneracy[codon[i]]:
                        basedex[b] = 1 / len(degeneracy[codon[i]])
                    else:
                        basedex[b] = 0
                pred.append(basedex)
            codonmix.append(pred)
        return list(zip(freq, codonmix))


if __name__ == "__main__":
    # the documentation is shocking
    # http://biopython.org/DIST/docs/api/Bio.SeqIO.AbiIO-module.html
    from Bio import SeqIO
    import json

    file = "example data/ACE-AA-088-01-55Â°C-BM3-A82_19C-T7-T7minus1.ab1"
    x = Trace.from_filename(file)
    print(x.QQC('CGT GAT TTT', '12NDT 6VHA 1TGG 1ATG'))
    # To Do...
    # calculate AA
    # give presents words like Tang and 22c
