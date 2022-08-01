#!/usr/bin/env python

###########################################################
###########################################################
###### Originally written by Habiballah Rahimi Eichi ######
###########################################################
###########################################################

import os
import sys
import argparse as ap
import logging
import re
from dateutil import tz
from datetime import datetime, date, timedelta
import pandas as pd
import gzip
from operator import itemgetter
from math import floor
from glob import glob
import pytz
import time
import subprocess as sp

logger = logging.getLogger(os.path.basename(__file__))
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_args():
    argparser = ap.ArgumentParser('freq Pipeline for GENEActiv')

    # Input and output parameters
    argparser.add_argument('--read-dir',
        help='Path to the input directory', required=True)
    argparser.add_argument('--output-dir',
        help='Path to the output directory', required=True)
    argparser.add_argument('--ext-mode')
    argparser.add_argument('--mtl-dir')
    return argparser

def main(args):
    # expand any ~/ in the directories
    read_dir = os.path.expanduser(args.read_dir)
    output_dir = os.path.expanduser(args.output_dir)
    ext_mode = os.path.expanduser(args.ext_mode)

    # perform sanity checks for inputs
    read_dir = check_input(read_dir)
    output_dir = check_output(output_dir)
    if read_dir is None or output_dir is None:
        return

    # logger output
    fh = logging.FileHandler(os.path.join(output_dir, 'freq.log'))
    logger.addHandler(fh)

    # run MATLAB
    run_matlab(read_dir, output_dir, args.mtl_dir, ext_mode)

# Run MATLAB
def run_matlab(read_dir, output_dir, mtl_dir,ext_mode):
    try:
        logger.info('Running matlab')
        matlab_path = "addpath('{matlab_dir}');".format(matlab_dir=mtl_dir)
        sub_cmd = "freq('{READ_DIR}','{OUTPUT_DIR}','{EXT_MODE}')".format(READ_DIR=read_dir,
            OUTPUT_DIR=output_dir, EXT_MODE=ext_mode)
        sub_cmd = wrap_matlab(sub_cmd)

        if mtl_dir:
            sub_cmd = matlab_path + sub_cmd

        cmd = ['matlab', '-nojvm', '-nodisplay', '-nosplash', '-r', sub_cmd]
        sp.check_call(cmd, stderr=sp.STDOUT)
    except Exception as e:
        logger.error(e)

def wrap_matlab(cmd):
    return 'try; {0}; catch; err = lasterror; disp(err.message); quit(1); end; quit();'.format(cmd)

# Exit program if the input directory does not exist.
def check_input(read_dir):
    if os.path.exists(read_dir):
        return os.path.join(read_dir, 'mtl1')
    else:
        logger.error('%s does not exist.' % read_dir)
        return None

# Exit program if the output directory does not exist.
def check_output(output_dir):
    if os.path.exists(output_dir):
        output_dir = os.path.join(output_dir, 'mtl2')
        if os.path.exists(output_dir):
#            clean_output_dir(output_dir)
            return output_dir
        else:
            try:
                os.mkdir(output_dir)
                return output_dir
            except Exception as e:
                logger.error('Could not create %s' % output_dir)
                return None
    else:
        logger.error('%s does not exist.' % output_dir)
        return None

def clean_output_dir(output_dir):
    logger.warn('Cleaning out output directory %s' % output_dir)
    file_name = 'mtl2_*'
    file_pattern = os.path.join(output_dir, file_name)
    for match in glob(file_pattern):
        logger.warn('Removing file %s' % match)
        os.remove(match)

if __name__ == '__main__':
    parser = parse_args()
    args = parser.parse_args()
    main(args)
