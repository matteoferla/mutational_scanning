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
if __name__ == "__main__":
    from Bio import SeqIO
    import json
    handle = open("example data/ACE-AA-088-01-55Â°C-BM3-A82_19C-T7-T7minus1.ab1", "rb")
    for record in SeqIO.parse(handle, "abi"):
        print(json.dumps(record.annotations, sort_keys=True, indent=4))