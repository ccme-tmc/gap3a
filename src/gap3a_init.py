#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Python version of gap_init.
'''

import os
from argparse import ArgumentParser, RawTextHelpFormatter#, RawDescriptionHelpFormatter


def gap_init():
    '''Initialize all necessary inputs for gap3a
    '''
    parser = ArgumentParser(description=__doc__, \
            formatter_class=RawTextHelpFormatter)
    
    parser.add_argument("-b", dest="batch_flag", action="store_true", \
        help="set batch mode, supress the interactive step")
    parser.add_argument('-c', dest='cmplx_flag', action="store_true", \
        help="wave functions are complex (no inversion symmetry present)")
    parser.add_argument("-t", dest="task", type=str, default="gw", \
        help="task name")
    parser.add_argument("-f", dest="case", type=str, default=".", \
        help="case name")
    parser.add_argument("-d", dest="gwdir", type=str, default="./gw", \
        help="directory to run the gw calculations")
    parser.add_argument("--sv", dest="symvec", type=int, default=1, \
        help="control how to use symmetry when generating k-mesh used for GW calculations")

    
    # initialize options as 'opts'
    opts = parser.parse_args()
    if opts.file == ".":
        case = os.path.basename(os.environ["PWD"])
    else:
        case = opts.case


if __name__ == "__main__":
    gap_init()
