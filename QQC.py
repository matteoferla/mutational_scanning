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
        '''
        Class method to return a Trace instance, but from a filename.
        >>> data=Trace.from_filename('sample.ab1')
        :param file: a file name
        :return: Trace instance
        '''
        handle = open(file, 'rb')
        record = SeqIO.read(handle, "abi")
        return cls(record)

    def find_peak(self, target_seq, strict=True):
        '''
        Find the index of the first base of a given sequence in the trace.
        :param target_seq: a string of bases.
        :param strict: boolean, dies if multiple hits.
        :return: the index of peak_id or peak_index of the first base.
        '''
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
        '''
        Same as find_peak, but give the index of the base after the given sequence.
        :param target_seq: a string of bases.
        :param strict: boolean, dies if multiple hits.
        :return: the index of peak_id or peak_index of the next (last+1) base along
        '''
        target_seq = re.sub('[^ATGC]', '', target_seq.upper())
        return self.find_peak(target_seq, strict)+len(target_seq)+1

    def get_intensities(self, index,wobble=0.20):
        '''
        Get intensities at given peak index.
        :param index: a peak index
        :param wobble: window to look in (default .2)
        :return: dictionary
        '''
        span = round(len(self.A)/len(self.peak_index))
        doubleindex=self.peak_index[index]
        return {base: max([getattr(self,base)[doubleindex+i] for i in range(-round(span*wobble), round(span*wobble))]) for base in self.bases}

    def QQC(self,location, *args, **kwargs):
        if isinstance(location, str):
            location = self.find_peak_after(location)
            ## go..
        QQC([self.get_intensities(location)],*args, **kwargs)


class QQC:
    def __init__(self,trace,location, scheme='NNK'):
        pass

    @staticmethod
    def from_trace(trace, location, *args, **kwargs):
        '''
        An alternate route. QQC instance comes from Trace(...).QQC()
        :param trace: Trace instance
        :param location: location of first peak of three
        :return: QQC instance
        '''
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
    print(x.peak_index)
    print(x.get_intensities(x.find_peak_after('CGT GAT TTT')))
    # To Do...
    # get base frequency at given peak
    # QQC class
    # given input CGT GAT TTT NNK
    # calculate ratios
    #
