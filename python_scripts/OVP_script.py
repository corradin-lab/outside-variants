#!/lab/corradin_data/FOR_AN/anaconda3/bin/python

"""
Authors: Parker Hall and An Hoang
"""

import sys
import itertools
import numpy as np
import pandas as pd
from collections import defaultdict, namedtuple
from functools import wraps

import argparse

# for file io
import os
import errno
from contextlib import contextmanager
import logging
from utils import cd, make_working_dir, remove_folder, make_hash
from pipelines import Pipeline

# TODO: turn this into a method of args object


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='outside var pipeline')
    parser.add_argument('input_folder_path',
                        help='absolute path of input folder')
    parser.add_argument(
        'case_gen', help='name chromosome-specific case .gen file, with path specified using argument PATH in init file. Example: ALL_MS_impute2_chr20.gen ')
    parser.add_argument(
        'control_gen', help='chromosome-specific control .gen file, with path specified using argument PATH in init_file. Example: ALL_controls_58C_NBS_WTC2_impute2_chr20.gen')
    parser.add_argument(
        'SNP_pairs', help='file with each line having GWAS_rsid and outside_rsid, separated by a delimiter. Example: SNPpairs_SAMPLE ')
    parser.add_argument(
        'init_file', help='file with init arguments. Formatted in the form: keyword,value separated by tab')
    parser.add_argument('output_folder', help='path for output folder for this run. WARNING: If --override flag is present, and there exists a folder with the same name at the path provided, the folder will be deleted before the pipeline runs. If unsure, do not supply --override.')
    parser.add_argument('unique_identifier', help='unique identifier for this run. Will be used for creating hash, and names of folders and files. This parameter is only important if you are using the orchestrator script. If you are calling super_pipeline script directly, you can make this the same as \'output_folder\'')
    parser.add_argument('--override', action='store_true',
                        help='if this argument is present, override any current folder with the same unique_identifier')
    parser.add_argument('--debug_mode', action='store_true', help="if this flag is present, set the logging engine to 'DEBUG' to print out more information during the run")

    parser.add_argument('--one_pair', "-op", nargs="?", default=None,
                        metavar="\{GWAS\}_\{Outside\}", help='supply SNP pair for directory structures (to be used with orchestrator.py only)')
    parser.add_argument('--iter', type=int, default=None,
                        help='if this argument is present, use the number of iterations from this argument instead of what specified in init file')

    parser.add_argument('--case', type=int, default=None, help='number of total cases. Usually will get this number from case sample file. Use for pipelines that do not specify sample files. WARNING: if both this flag and the sample case is present, it will cause an error and terminate the program.')

    parser.add_argument('--control', type=int, default=None, help='number of total control. Usually will get this number from case sample file. Use for pipelines that do not specify sample files. WARNING: if both this flag and the sample case is present, this will cause an error and terminate the program.')

    args = parser.parse_args()
    if args.debug_mode:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logging.info("Starting outside variant pipeline analysis")
    file1 = args.case_gen
    file2 = args.control_gen
    pairing = args.SNP_pairs
    init_file = args.init_file
    p_file = args.output_folder
    override_folder = args.override


    odds_file = ""

    logging.info("Initializing pipeline. This might take a few seconds.")
    args.exec_dir = os.getcwd()
    with cd(args.input_folder_path):
        pipe = Pipeline.init_from_file(
            init_file, file1, file2, pairing, p_file, odds_file, args)
    logging.info("Making output directory")
    working_dir = make_working_dir(p_file, override_folder)
    pipe.working_dir = working_dir
    pipe.p_value_filename = p_file.split("/")[-1]
    pipe.hash = make_hash(args.input_folder_path, init_file,
                          file1, file2, pairing, args.unique_identifier)

    with cd(args.input_folder_path):
        pipe.read_input_files()
    logging.info("Running pipeline...")
    with cd(pipe.working_dir):
        pipe.run()
