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
        ## sanitise
        targetdex = re.findall(target_seq, self.peak_id)
        if not targetdex:
            raise ValueError('target_seq not found!')
        elif len(targetdex) > 1:
            if strict:
                raise ValueError('Ambiguous sequence given!')
            else:
                warn('Ambiguous sequence given!')
        i = self.peak_id.find(target_seq) + len(target_seq)
        ## analyse


if __name__ == "__main__":
    # the documentation is shocking
    # http://biopython.org/DIST/docs/api/Bio.SeqIO.AbiIO-module.html
    from Bio import SeqIO
    import json

    file = "example data/ACE-AA-088-01-55Â°C-BM3-A82_19C-T7-T7minus1.ab1"
    x = Trace.from_filename(file)
    x.find_peak('ACGTGATTTT')
