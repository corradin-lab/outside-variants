import logging
import os

from MTC import MtcTable, StepwiseFilter
from collections import defaultdict, namedtuple

import pandas as pd
import numpy as np
from time import time, sleep
# package to saving and loading dict
import pickle
from tqdm.auto import tqdm
from utils import atomic_write, timeit, open_file_or_string, possible_genotypes, randomizer, cd, swap_attr

logging = logging.getLogger('OVP')

def process_cli_args(pipe, cli_args):
    # changing pipeline attributes based on command line arguments

    pipe.one_pair = cli_args.one_pair

    if pipe.one_pair:
        GWAS_rsID, *outside_rsID = pipe.one_pair.split("_")

        #for outside rsIDs in the form c6_pos1234:
        outside_rsID = "_".join(outside_rsID)

        pipe.rsid_pair_dir = os.path.join(
            cli_args.exec_dir, "{}/{}".format(GWAS_rsID, outside_rsID))

    abs_input_folder_path = os.path.abspath(cli_args.input_folder_path)

    # cli_args.input_folder_path
    pipe.args_input_folder_path = abs_input_folder_path
    pipe.input_folder_path = ("/").join(abs_input_folder_path.split("/")[:-1])
    logging.debug(f"absolute input folder path: {pipe.input_folder_path}")
    if cli_args.iter:
        pipe.iterations = cli_args.iter

    # for case/control, if one of them is missing while the other is supplied.
    # raise an error
    if (not cli_args.case) == (not cli_args.control):
        pipe.total_case_unfiltered = cli_args.case
        pipe.total_control_unfiltered = cli_args.control
    else:
        raise ValueError(
            "case and control must be both specified or both not specified")
    return pipe

class Pipeline:
    """
    ABSTRACT CLASS: initializing an instance of Pipeline will not work
                                                                    instead, initialize one of the classes below that inherit from Pipeline
    """

    def __init__(self, pair, rand, p, gs, itr, typ, rsid, snp, delim, or_calc, skip, path, partial, check, covs, high, low, NA):
        self.pair_filename = pair
        self.rand_filename = rand
        self.p_value_filename = p
        self.col_cut = gs
        self.iterations = itr
        self.is_triples = typ
        self.rsid_col = rsid
        self.snp_col = snp
        self.delim = delim
        self.odds_ratio_type = or_calc
        self.skip = "<DEL>" if skip is None else skip
        self.filepath = path
        self.output_partial = partial
        self.check_start = check

        self.high_lim = .9 if high is None else high
        self.low_lim = .3 if low is None else low

        self.NA = "NA" if NA is None else NA
        self.sorted_NA = "".join(sorted(self.NA))

        self.odds_ratio_index = 6
        self.low_p_index = 7
        self.high_p_index = 8
        self.total_index = 9
        self.int_indices = set([0, 1, 3, 4, 9])

        self.num_decimals = 4

        self.output_odds = bool(self.rand_filename)

        self.covariates = self.read_covariate_file(covs)

        # TO DO: working
        #

        if self.odds_ratio_type == 3:
            import sklearn.linear_model as md
            self.model = md.LogisticRegression(solver='liblinear')

        self.RANDOMIZER_CLEAN_DICT = {}  # cache of cleaned strings for faster lookup

        if self.output_partial:
            # create output files to later be appended to
            if self.output_odds:
                f = open(self.rand_filename, "w")
                f.close()

    @staticmethod
    def init_from_file(init_file, file1, file2, pair, p, rand, cli_args):
        """
        initializes pipeline based off of init_file
        """
        MAIN_PIPELINE_KEYWORD_DICT = {
            "PIPE_TYPE": 0,
            "GS": 1,
            "ITER": 2,
            "TRIPS": 3,
            "RSID": 4,
            "SNP": 5,
            "DELIM": 6,
            "OR_CALC": 7,
            "SKIP": 8,
            "PATH": 9,
            "PARTIAL": 10,
            "CHECK": 11,
            "SAMPLE": 12,
            "COVS": 13,
            "CUTOFF": 14,
            "NA": 15,
        }

        SAMPLE_FILE_KEYWORD_DICT = {
            "LINES_TO_SKIP": 0,
            "STATUS_COL_NAME":1,
            "CASE_INT":2,
            "CONTROL_INT":3,
        }

        DELIM_KEYWORD_DICT = {
            "DELIM_GEN":0,
            "DELIM_PAIRING": 1,
            "DELIM_SAMPLE":2
        }

        GET_MTC_KEYWORD_DICT = {
            "MTC_FILE_PATH": 0,
            "MTC_THRESHOLD_COL_LABEL": 1,
            "Z_THRESHOLD": 2
        }

        MAKE_MTC_KEYWORD_DICT = {
            "MAKE_MTC": 0
        }

        INT_KEYWORDS = set(["GS", "ITER", "TRIPS", "RSID", "SNP", "OR_CALC",
                            "PARTIAL", "CHECK", "Z_THRESHOLD", "MAKE_MTC"])
        DELIM_KEYWORDS = DELIM_KEYWORD_DICT.keys()

        def keywords_to_lists(*keywords_dicts):

            # a list of values for each of the keyword dict
            init_lists = [[None] * len(d) for d in keywords_dicts]
            logging.info(f"Reading init file at {os.getcwd()}/{init_file}")
            try:
                f = open(init_file, "r")
            except FileNotFoundError:
                raise FileNotFoundError(
                    "Cannot read init file, please check if they exist in the directory {} and in the right format. Also, make sure they don't have hidden file extensions like .txt".format(os.getcwd()))

            for raw_line in f:

                clean_line = raw_line.strip()

                if clean_line:
                    keyword, value = clean_line.split("\t")
                    keyword = keyword.upper()

                    if keyword in INT_KEYWORDS:
                        value = int(value)
                    elif keyword in DELIM_KEYWORDS :
                        if value == "TAB":
                            value = "\t"
                        elif value == "SPACE":
                            value = " "
                        else:
                            raise ValueError(
                                "Delimiter type not supported. Please try 'SPACE' or 'TAB'.")

                    # find the keyword in each of the keyword dicts
                    try:
                        dict_index, found_dict = next(
                            (i, d) for i, d in enumerate(keywords_dicts) if keyword in d)

                        keyword_index = found_dict[keyword]
                        init_lists[dict_index][keyword_index] = value
                    # raise error if cannot find the argument in any of the
                    # dicts passed in
                    except StopIteration:
                        raise KeyError(
                            "Argument keyword {} not recognized. Check example init file.".format(keyword))

            return init_lists

        def get_mtc_with_error_logging(pipe, get_mtc_init_list):
            logging.debug("GETTING MTC TABLE")
            if all(get_mtc_init_list):
                try:
                    pipe.get_mtc_table_and_stepwise_filter(*get_mtc_init_list)
                    pipe.run = pipe.run_with_permutation()
                except ValueError as e:
                    raise e("Failed to get mtc table and stepwise filter object. Check that all the following arguments in init file are correct: {}".format(
                        GET_MTC_KEYWORD_DICT.keys()))
            else:
                raise ValueError("MISSING ARGUMENT IN INIT FILE. Failed to get MTC table. Make sure the following arguments in init file are not empty: {}".format(
                    GET_MTC_KEYWORD_DICT.keys()))

        def make_mtc_with_error_logging(pipe, make_mtc_init_list, case_sample_file, control_sample_file, case_gen, control_gen):
            """modifying the run function of the pipeline so it can make MTC table or run permutation

            Args:
                    pipe (TYPE): Description
                    make_mtc_init_list (TYPE): Description
            """
            logging.debug("CHECKING SAMPLE FILES")
            # if both sample file and --case/--control flags are specified,
            # raise an error

            if any([case_sample_file, control_sample_file]) and any([pipe.total_control_unfiltered]):
                raise ValueError(
                    " The specified --case/--control flags conflict with the sample files, please choose one or the other and not both")

            logging.debug("MAKING MTC TABLE")

            if all([case_sample_file, control_sample_file]):
                try:
                    # num = length of file minus header
                    pipe.total_case_unfiltered = pd.read_table(
                        case_sample_file).shape[0] - 1
                    pipe.total_control_unfiltered = pd.read_table(
                        control_sample_file).shape[0] - 1
                except FileNotFoundError:
                    raise FileNotFoundError(
                        "Cannot read sample files, please check if they exist in the directory {} and in the right format. Also, make sure they don't have hidden file extensions like .txt".format(pipe.input_folder_path))

            if all([pipe.total_case_unfiltered, pipe.total_control_unfiltered]):
                pipe.run = pipe.run_create_mtc_table(
                    pipe.total_control_unfiltered, pipe.total_control_unfiltered, case_gen, control_gen)
            else:
                raise ValueError(
                    "Cannot find total case and total control, please supply either sample files in init file with the keyword SAMPLE, or in the command line using \'--case\' flag and \'--control\' flag")

        # the order that you pass the dictionaries in matters. It determines
        # the order of the lists you get out
        main_pipeline_init_list, get_mtc_init_list, make_mtc_init_list, delim_init_list = keywords_to_lists(
            MAIN_PIPELINE_KEYWORD_DICT, GET_MTC_KEYWORD_DICT, MAKE_MTC_KEYWORD_DICT, DELIM_KEYWORD_DICT)

        pipe_type, gs, itr, typ, rsid, snp, delim, or_calc, skip, path, partial, check, sample, covs, cuts, NA = main_pipeline_init_list

        if sample != None:
            case, control = sample.split(",")
        else:
            case, control = None, None

        if cuts != None:
            high, low = cuts.split(",")
            high = float(high)
            low = float(low)
        else:
            high, low = None, None

        # create correct type of pipeline
        if pipe_type == "CC":
            if or_calc == 3 and sample == None:
                raise ValueError(
                    "Cannot calculate logistic odds ratio without sample files")
            pipe = Case_Control_Pipeline(file1, file2, pair, rand, p, gs, itr, typ, rsid, snp,
                                         delim, or_calc, skip, path, partial, check, case, control, covs, high, low, NA)

        elif pipe_type == "COMB":
            pipe = Combined_Pipeline(file1, file2, pair, rand, p, gs, itr, typ, rsid,
                                     snp, delim, or_calc, skip, path, partial, check, covs, high, low, NA)
        elif pipe_type == "TRANS_CC":
            if or_calc == 3 and sample == None:
                raise ValueError(
                    "Cannot calculate logistic odds ratio without sample files")
            pipe = Trans_Case_Control_Pipeline(file1, file2, pair, rand, p, gs, itr, typ, rsid,
                                               snp, delim, or_calc, skip, path, partial, check, case, control, covs, high, low, NA)

        elif pipe_type == "TRANS_COMB":
            pipe = Trans_Combined_Pipeline(file1, file2, pair, rand, p, gs, itr, typ,
                                           rsid, snp, delim, or_calc, skip, path, partial, check, covs, high, low, NA)
        else:
            raise ValueError(
                "Pipeline type not supported. Try 'CC', 'COMB', 'TRANS_CC', or 'TRANS_COMB'.")

        logging.debug("before processing cli args")

        pipe = process_cli_args(pipe, cli_args)

        # store in global variable for stepwise filter
        global PIPE_TYPE
        PIPE_TYPE = pipe.__class__

        # make sure that you cannot get contradictory arguments of making and
        # getting MTC table at the same time
        if any(get_mtc_init_list) and any(make_mtc_init_list):
            raise ValueError("Cannot have arguments to get MTC table (arguments: {}) and make MTC table (arguments: {}) at the same time ".format(
                GET_MTC_KEYWORD_DICT.keys(), MAKE_MTC_KEYWORD_DICT.keys()))

        # has to either make mtc file or use a ready-made mtc file
        if all(get_mtc_init_list):
            pipe.make_mtc = False
            get_mtc_with_error_logging(pipe, get_mtc_init_list)
        elif all(make_mtc_init_list):
            pipe.make_mtc = True
            make_mtc_with_error_logging(
                pipe, make_mtc_init_list, case, control, file1, file2)
        else:
            raise ValueError(f"Need to either specify arguments: {list(GET_MTC_KEYWORD_DICT.keys())} or {list(MAKE_MTC_KEYWORD_DICT.keys())} in init file {init_file}")

        #process delimiter args
        pipe.gen_delim, pipe.pairing_delim, pipe.sample_delim = delim_init_list
        return pipe

    def get_mtc_table_and_stepwise_filter(self, mtc_file_path, mtc_threshold_col_label, z_threshold):

        mtc_file_exists = os.path.isfile(mtc_file_path)
        # TODO: FINISH ELSE CLAUSE
        if mtc_file_exists:
            self.mtc_table = MtcTable(mtc_file_path, mtc_threshold_col_label)
            self.stepwise_filter = StepwiseFilter(self.mtc_table, z_threshold)
        else:
            raise ValueError(
                "Cannot find MTC table at designated path: " + mtc_file_path + f".Currently in path {os.getcwd()}")
        # else:
        #   self.mtc_table = Mtc_table

    def triplicate_converter(self, allele_list, pos1, pos2):
        """
        takes in allele_list in gen triplicate format and
        return new_list in letter format
        """
        all_geno = possible_genotypes(pos1, pos2)
        new_list = []
        for i in range(0, len(allele_list), 3):
            geno_1 = float(allele_list[i])
            geno_2 = float(allele_list[i + 1])
            geno_3 = float(allele_list[i + 2])
            if geno_1 >= self.high_lim and geno_2 < self.low_lim and geno_3 < self.low_lim:
                new = all_geno[0]
            elif geno_2 >= self.high_lim and geno_1 < self.low_lim and geno_3 < self.low_lim:
                new = all_geno[1]
            elif geno_3 >= self.high_lim and geno_1 < self.low_lim and geno_2 < self.low_lim:
                new = all_geno[2]
            else:
                new = self.NA
            new_list.append(new)
        return new_list

    def clean_item_cached(self, item):
        """
        cleans item string. caches results for constant lookups of known items
        """
        return self.RANDOMIZER_CLEAN_DICT.get(item, None) or self.clean_item(item)

    def clean_item(self, item):
        orig_item = item
        item = "".join(sorted(item))
        if item == self.sorted_NA:
            item = self.NA
        if " " in item:
            item = item.replace(" ", "")
        self.RANDOMIZER_CLEAN_DICT[orig_item] = item
        return item

    def read_covariate_file(self, file):
        """
        reads single-column file, adding all elements to list
        returns list
        """
        covs = []
        if file != None:
            f = open(file, "r")
            for line in f:
                line = line.strip().lower()
                covs.append(line)
            f.close()
        return covs

    def extract_cov_info(self, line, filename):
        cov_indices = []
        cov_not_found = []
        for covariate in self.covariates:
            cov_index = self.index_getter(line, covariate)
            if cov_index > -1:
                cov_indices.append(cov_index)
            else:
                cov_not_found.append(covariate)

        if len(self.covariates) > 0:
            if len(cov_not_found) == 0:
                logging.info("ALL COVARIATES FOUND")
            elif len(cov_indices) == 0:
                raise LookupError(
                    "WARNING: none of the specified covariates were found in " + filename)
            else:
                logging.warning("WARNING: the following covariates were not found in " +
                                filename + ": " + str(cov_not_found))

        return cov_indices

    def column_getter(self, line, col_num):
        """
        returns item in the 'col_num'th column in tab-delimited line (1 indexed)
        """

        copy_line = line.strip()
        for i in range(col_num - 1):
            tab_index = copy_line.index(self.delim)
            copy_line = copy_line[tab_index + 1:]
        if self.delim in copy_line:
            #logging.debug("found delim")
            final_index = copy_line.index(self.delim)
            out = copy_line[:final_index]
        else:
            #logging.debug("did not find delim")
            out = copy_line
        return out

    def column_cutter(self, line, col_num):
        """
        cuts tab-delimited line before col_num, returns pieces as a tuple
        (first piece includes trailing tab)
        """
        after = line.strip()
        for i in range(col_num - 1):
            tab_index = after.index(self.delim)
            after = after[tab_index + 1:]
        ind = line.find(after)
        before = line[:ind]
        return (before, after)

    def index_getter(self, string, target):
        """
        takes in string delimited by delimiter
        returns index of target string in list of strings split by delimiter
        (zero-indexed)
        """
        string = string.strip()
        line_list = string.split(self.delim)
        for i in range(len(line_list)):
            if line_list[i].lower() == target:
                return i
        return -1

    def check_format_single_file(self, filename):
        """
        checks that correct genotype start column has been specified
        """
        check = open(filename, "r")

        line_counter = 1
        for line in check:
            dashes = self.column_getter(line, self.col_cut - 1)
            SNP = self.column_getter(line, self.col_cut)

            bad_dashes = dashes.isalpha()
            bad_SNP = not SNP.isalpha()
            if self.is_triples:
                bad_dashes = not bad_dashes and dashes != self.skip
                bad_SNP = not bad_SNP

            if bad_dashes or bad_SNP:
                check.close()
                raise ValueError(
                    "Check " + filename + ". Incorrect SNP start column on line " + str(line_counter) + ".")
            line_counter += 1

        check.close()

    @timeit
    def pairing_maker(self):
        """
        reads pairing file and creates GWAS=>outside rsID dictionary,
        outside=>GWAS rsID dictionary, GWAS rsID set, and outside rsID set
        returns these four data structures
        """

        self.G2O = defaultdict(list)
        self.O2G = defaultdict(list)
        self.GWAS_set = set()
        self.outside_set = set()
        self.all_SNP_pairs_in_pairing_file_set = set()

        with open_file_or_string(self.pair_filename, "r") as pairing_file:
            for line in pairing_file:
                GWAS_rsID = self.column_getter(line, 1)
                outside_rsID = self.column_getter(line, 2)
                self.GWAS_set.add(GWAS_rsID)
                self.outside_set.add(outside_rsID)

                # record the original SNP pairs in the pairing file
                self.all_SNP_pairs_in_pairing_file_set.add(
                    (GWAS_rsID, outside_rsID))

                self.G2O[GWAS_rsID].append(outside_rsID)
                self.O2G[outside_rsID].append(GWAS_rsID)

        if len(self.G2O) == 0 or len(self.O2G) == 0:
            raise IOError("WARNING: SNP pair file is EMPTY")

        return self.G2O, self.O2G, self.GWAS_set, self.outside_set

    def odds_ratio_calculator(self, case, control, case_total, control_total):
        """
        calculates odds ratio using formula specified by type of pipeline
        """
        try:
            if self.odds_ratio_type == 1:
                case_odds = case / (case_total - case)
                control_odds = control / (control_total - control)
            else:
                case_odds = case / case_total
                control_odds = control / control_total

            odds_ratio = case_odds / control_odds
            if odds_ratio == 0:
                odds_ratio = "NA"

            return odds_ratio
        except ZeroDivisionError:
            return "NA"

    @timeit
    def log_odds_calculator(self, case, control):
        num_cases = case["total"]
        num_controls = control["total"]

        case_NA = case["NA_indices"]
        control_NA = list(map(lambda x: x + num_cases, control["NA_indices"]))

        NA_indices = case_NA + control_NA

        labels = np.zeros(num_cases + num_controls)
        np.put(labels, range(num_cases), [1])
        labels = np.delete(labels, NA_indices)

        covs_no_NA = np.delete(self.covariate_data, NA_indices, 0)

        odds_ratio_dict = {}

        non_genotypes = set(["total", "NA_indices", "alleles"])

        for genotype in case:
            if genotype not in non_genotypes:
                data = np.zeros(num_cases + num_controls)

                case_hits = case[genotype]
                control_hits = list(
                    map(lambda x: x + num_cases, control[genotype]))

                if len(case_hits) == 0 or len(control_hits) == 0:
                    odds_ratio = "NA"
                else:
                    np.put(data, case_hits, [1])
                    np.put(data, control_hits, [1])

                    data = np.delete(data, NA_indices).reshape((-1, 1))
                    data_with_covs = np.hstack((data, covs_no_NA))

                    # logging.debug("Features: ")
                    # logging.debug(data_with_covs)

                    # logging.debug("labels:")
                    # logging.debug(labels)

                    # logging.debug("----")
                    #model = self.model.LogisticRegression().fit(data_with_covs, labels)

                    model = self.model.fit(data_with_covs, labels)

                    beta = model.coef_[0][0]
                    odds_ratio = np.exp(beta)

                odds_ratio_dict[genotype] = odds_ratio
        return odds_ratio_dict

    @timeit
    def GWAS_OR_calc(self):
        """
        calculates GWAS odds ratio and creates GWAS rsID => odds ratio dictionary
        returns this dictionary
        """
        self.GWAS_OR_dict = {}

        for rsid in self.case_GWAS:
            if rsid in self.control_GWAS:
                case_inner_dict = self.case_GWAS[rsid]
                control_inner_dict = self.control_GWAS[rsid]

                genotypes = case_inner_dict["alleles"]

                # compute logistic regression for odds ratio
                if self.odds_ratio_type == 3:
                    odds_ratio_dict = self.log_odds_calculator(
                        case_inner_dict, control_inner_dict)
                    #logging.debug("odds_ratio_dict: ", odds_ratio_dict)
                    self.GWAS_OR_dict[rsid] = odds_ratio_dict
                    # logging.debug(self.GWAS_OR_dict)
                # extract information for linear odds ratio
                elif self.odds_ratio_type == 2:
                    AA, AB, BB = sorted(genotypes)

                    num_AA_cases = float(len(case_inner_dict[AA]))
                    num_AB_cases = float(len(case_inner_dict[AB]))
                    num_BB_cases = float(len(case_inner_dict[BB]))

                    num_AA_controls = float(len(control_inner_dict[AA]))
                    num_AB_controls = float(len(control_inner_dict[AB]))
                    num_BB_controls = float(len(control_inner_dict[BB]))

                    num_A_cases = num_AA_cases * 2 + num_AB_cases
                    num_B_cases = num_BB_cases * 2 + num_AB_cases

                    num_A_controls = num_AA_controls * 2 + num_AB_controls
                    num_B_controls = num_BB_controls * 2 + num_AB_controls

                    if num_A_cases + num_A_controls <= num_B_cases + num_B_controls:
                        first_cases, first_controls = num_A_cases, num_A_controls
                        second_cases, second_controls = num_B_cases, num_B_controls
                        order = AB[0] + "/" + AB[1]
                    else:
                        first_cases, first_controls = num_B_cases, num_B_controls
                        second_cases, second_controls = num_A_cases, num_A_controls
                        order = AB[1] + "/" + AB[0]

                    odds_ratio = self.odds_ratio_calculator(
                        first_cases, second_cases, first_controls, second_controls)
                    self.GWAS_OR_dict[rsid] = (odds_ratio, order)

                # extract information for nonlinear odds ratio
                else:
                    case_total = float(
                        case_inner_dict["total"] - len(case_inner_dict["NA_indices"]))
                    control_total = float(
                        control_inner_dict["total"] - len(control_inner_dict["NA_indices"]))

                    new_map = {}
                    for geno in genotypes:
                        case_list = case_inner_dict[geno]
                        control_list = control_inner_dict[geno]

                        num_cases = len(case_list)
                        num_controls = len(control_list)
                        odds_ratio = self.odds_ratio_calculator(
                            num_cases, num_controls, case_total, control_total)
                        new_map[geno] = odds_ratio
                    self.GWAS_OR_dict[rsid] = new_map

        return self.GWAS_OR_dict

    @timeit
    def combine_dicts(self):
        """
        combines GWAS and outside genetic info for case and control groups
        """
        self.case_combined = self.combine_dicts_single_pair(
            self.case_GWAS, self.case_outside)
        self.control_combined = self.combine_dicts_single_pair(
            self.control_GWAS, self.control_outside)

        return self.case_combined, self.control_combined

    def combine_dicts_single_pair(self, GWAS, outside):
        """
        combines GWAS and outside genetic info
        returns (GWAS rsID, outside rsID) => (GWAS allele, outside allele) => count dictionary
        """
        combined = {}
        for GWAS_key in self.G2O:
            if GWAS_key in GWAS:
                small_GWAS_dict = GWAS[GWAS_key]
                GWAS_NA_indices = small_GWAS_dict["NA_indices"]
                GWAS_total = small_GWAS_dict["total"]
                GWAS_alleles = small_GWAS_dict["alleles"]
                for outside_key in self.G2O[GWAS_key]:
                    if outside_key in outside:
                        outside_list, outside_alleles = outside[outside_key]

                        new_map = defaultdict(list)
                        new_map["NA_indices"] = GWAS_NA_indices[:]
                        new_map["total"] = GWAS_total

                        for GWAS_allele in GWAS_alleles:
                            GWAS_indices = small_GWAS_dict[GWAS_allele]
                            for GWAS_index in GWAS_indices:
                                outside_allele = outside_list[GWAS_index]
                                if outside_allele == self.NA:
                                    new_map["NA_indices"].append(GWAS_index)
                                else:
                                    new_key = (GWAS_allele, outside_allele)
                                    new_map[new_key].append(GWAS_index)

                            for poss_outside in outside_alleles:
                                check = (GWAS_allele, poss_outside)
                                new_map.setdefault(check, [])

                        combined[(GWAS_key, outside_key)] = new_map
        return combined

    @timeit
    def nonrand_dict_maker(self):
        """
        creates real GWAS-outside dictionary containing case/control and odds ratio info
        """

        # added to save dictionaries to load instead of recreating every time
        # odds_ratio_dict_hashed_filename = f"odds_ratio_dict_{unique_hash}"
        odds_ratio_dict_hashed_filename = "odds_ratio_dict_{}".format(
            self.hash)
        # p_value_dict_hashed_filename = f"p_value_dict_{unique_hash}"
        p_value_dict_hashed_filename = "p_value_dict_{}".format(self.hash)

        try:
            self.p_value_dict = self.load_dict(p_value_dict_hashed_filename)
            logging.debug("Loaded nonrand dict from file")
        except:
            logging.debug(
                "Could not find saved dict file, creating dicts from scratch")
            nonrandom_map = {}
            if self.output_odds:
                self.odds_ratio_dict = {}

            for G_O_rsid in self.case_combined:
                if G_O_rsid in self.control_combined:
                    geno_to_case_count = self.case_combined[G_O_rsid]
                    geno_to_control_count = self.control_combined[G_O_rsid]

                    if self.odds_ratio_type == 3:
                        odds_ratio_dict = self.log_odds_calculator(
                            geno_to_case_count, geno_to_control_count)
                        # logging.debug(odds_ratio_dict)
                        for genotype in odds_ratio_dict:
                            odds_ratio = odds_ratio_dict[genotype]
                            case_count = len(geno_to_case_count[genotype])
                            control_count = len(
                                geno_to_control_count[genotype])

                            case_total = geno_to_case_count[
                                "total"] - len(geno_to_case_count["NA_indices"])
                            control_total = geno_to_control_count[
                                "total"] - len(geno_to_control_count["NA_indices"])

                            new_key = (G_O_rsid, genotype)
                            # logging.debug(self.GWAS_OR_dict)
                            nonrandom_map[new_key] = [case_count, case_total, "N/A",
                                                                              control_count, control_total, "N/A",
                                                                              odds_ratio, 0, 0, 0, self.GWAS_OR_dict[G_O_rsid[0]][genotype[0]]]
                            if odds_ratio == "NA":
                                for ind in [self.low_p_index, self.high_p_index, self.total_index]:
                                    nonrandom_map[new_key][ind] = "NA"

                            if self.output_odds:
                                self.odds_ratio_dict[new_key] = [odds_ratio]

                    # extract information for linear odds ratio
                    elif self.odds_ratio_type == 2:
                        del geno_to_case_count["total"]
                        del geno_to_control_count["total"]

                        del geno_to_case_count["NA_indices"]
                        del geno_to_control_count["NA_indices"]

                        GWAS_dict_cases = {}
                        outside_genos = set()
                        for key in geno_to_case_count:
                            GWAS, outside = key
                            outside_genos.add(outside)
                            if GWAS not in GWAS_dict_cases:
                                GWAS_dict_cases[GWAS] = {}
                            GWAS_dict_cases[GWAS][outside] = float(
                                len(geno_to_case_count[key]))

                        GWAS_dict_controls = {}
                        for key in geno_to_control_count:
                            GWAS, outside = key
                            outside_genos.add(outside)
                            if GWAS not in GWAS_dict_controls:
                                GWAS_dict_controls[GWAS] = {}
                            GWAS_dict_controls[GWAS][outside] = float(
                                len(geno_to_control_count[key]))

                        GWAS_OR, GWAS_order = self.GWAS_OR_dict[G_O_rsid[0]]

                        AA, AB, BB = sorted(outside_genos)

                        for GWAS_geno in GWAS_dict_cases:
                            inner_case = GWAS_dict_cases[GWAS_geno]
                            inner_control = GWAS_dict_controls[GWAS_geno]

                            num_A_cases = inner_case[AA] * 2 + inner_case[AB]
                            num_A_controls = inner_control[
                                AA] * 2 + inner_control[AB]

                            num_B_cases = inner_case[BB] * 2 + inner_case[AB]
                            num_B_controls = inner_control[
                                BB] * 2 + inner_control[AB]

                            if num_A_cases + num_A_controls <= num_B_cases + num_B_controls:
                                first_cases, first_controls = num_A_cases, num_A_controls
                                second_cases, second_controls = num_B_cases, num_B_controls
                                order = AB[0] + "/" + AB[1]
                            else:
                                first_cases, first_controls = num_B_cases, num_B_controls
                                second_cases, second_controls = num_A_cases, num_A_controls
                                order = AB[1] + "/" + AB[0]

                            try:
                                first_odds = first_cases / first_controls
                            except ZeroDivisionError:
                                first_odds = "NA"

                            try:
                                second_odds = second_cases / second_controls
                            except ZeroDivisionError:
                                second_odds = "NA"

                            odds_ratio = self.odds_ratio_calculator(
                                first_cases, second_cases, first_controls, second_controls)
                            new_key = (G_O_rsid, (GWAS_geno,))
                            nonrandom_map[new_key] = [first_cases, first_controls, first_odds,
                                                      second_cases, second_controls, second_odds,
                                                      odds_ratio, 0, 0, 0, order, GWAS_OR, GWAS_order]

                            if odds_ratio == "NA":
                                for ind in [self.low_p_index, self.high_p_index, self.total_index]:
                                    nonrandom_map[new_key][ind] = "NA"

                            if self.output_odds:
                                self.odds_ratio_dict[new_key] = [odds_ratio]

                    # extract information for nonlinear odds ratio
                    else:
                        case_total = float(geno_to_case_count[
                            "total"] - len(geno_to_case_count["NA_indices"]))
                        control_total = float(geno_to_control_count[
                            "total"] - len(geno_to_control_count["NA_indices"]))

                        del geno_to_case_count["total"]
                        del geno_to_control_count["total"]

                        del geno_to_case_count["NA_indices"]
                        del geno_to_control_count["NA_indices"]

                        for genotype in geno_to_case_count:
                            case_count = len(geno_to_case_count[genotype])
                            control_count = len(
                                geno_to_control_count[genotype])

                            if self.odds_ratio_type == 0:
                                case_odds = case_count / case_total
                                control_odds = control_count / control_total
                            else:
                                try:
                                    case_odds = case_count / \
                                        (case_total - case_count)
                                except ZeroDivisionError:
                                    case_odds = "NA"

                                try:
                                    control_odds = control_count / \
                                        (control_total - control_count)
                                except ZeroDivisionError:
                                    control_odds = "NA"

                            odds_ratio = self.odds_ratio_calculator(
                                case_count, control_count, case_total, control_total)

                            new_key = (G_O_rsid, genotype)
                            nonrandom_map[new_key] = [case_count, case_total, case_odds,
                                                      control_count, control_total, control_odds,
                                                      odds_ratio, 0, 0, 0, self.GWAS_OR_dict[G_O_rsid[0]][genotype[0]]]
                            if odds_ratio == "NA":
                                for ind in [self.low_p_index, self.high_p_index, self.total_index]:
                                    nonrandom_map[new_key][ind] = "NA"

                            if self.output_odds:
                                self.odds_ratio_dict[new_key] = [odds_ratio]

            self.p_value_dict = nonrandom_map
            self.save_dict(self.p_value_dict, self.working_dir,
                           p_value_dict_hashed_filename)
            # testing
            #self.save_dict(self.p_value_dict, self.working_dir, "p_value_dict")
            #self.save_dict(self.odds_ratio_dict, odds_ratio_dict_hashed_filename)

        return self.p_value_dict

# functions for saving dictionaries into pickle files for easy i/o--------
    def save_dict(self, dict_to_save, target_dir, filename):
        # for one pair (using the orchestrator), we use the common rsid pair
        # dir to save the dicts for multiple filtering stages.

        if self.one_pair:
            save_dir = self.rsid_pair_dir
        else:
            # for multiple pairs, (not using the orchestrator), we use the
            # working directory
            save_dir = os.getcwd()

        dict_filename = "objs/" + filename + ".pickle"
        with cd(save_dir):
            self.try_make_folder("objs")

            with atomic_write(dict_filename, binary=True) as f:
                pickle.dump(dict_to_save, f, protocol=pickle.HIGHEST_PROTOCOL)

            logging.debug("SAVED FILE {} IN DIRECTORY {}".format(
                dict_filename, save_dir))

    def try_make_folder(self, folder_name_string):
        try:
            # Create target Directory
            os.makedirs(folder_name_string + "/")
            logging.debug("Creating the {} folder".format(folder_name_string))
        except:
            logging.debug("{} folder already exist".format(folder_name_string))

    def load_dict(self, filename, load_dir=None):
        # for one pair (using the orchestrator), we use the common rsid pair
        # dir to save the dicts for multiple filtering stages.
        if self.one_pair:
            logging.debug("rsid_pair_dir: ", self.rsid_pair_dir)
            load_dir_default = self.rsid_pair_dir#self.rsid_pair_dir
        else:
            # for multiple pairs, (not using the orchestrator), we use the
            # working directory
            load_dir_default = os.getcwd()

        if load_dir:
            load_dir_default = load_dir

        dict_filename = "objs/" + filename + ".pickle"
        with cd(load_dir_default):
            logging.debug("looking for file: {} in directory: {}".format(
                dict_filename, load_dir))
            with open(dict_filename, "rb") as f:
                dict_to_load = pickle.load(f)
            logging.debug("loaded file: {} from directory: {}".format(
                dict_filename, load_dir))
        return dict_to_load


#-----------------------------------------------------------------
#

    @timeit
    def perform_randomizations_an(self):
        """
        randomizes the outside genetic data the specified number of times
        calculates new odds ratios for every randomization
        updates p value and randomization data structures
        """

        # needs: self.control_outside[outside_rsid], self.O2G[outside_rsid]

        def outside_rsid_processor(self, outside_rsid):
            if outside_rsid in self.control_outside:
                case_alleles, poss_geno = self.case_outside[outside_rsid]
                control_alleles = self.control_outside[outside_rsid][0]

                AA, AB, BB = sorted(poss_geno)
                for GWAS_rsid in self.O2G[outside_rsid]:
                    GWAS_rsid_processor(
                        self, GWAS_rsid, case_alleles, poss_geno, control_alleles, AA, AB, BB)

        def GWAS_rsid_processor(self, GWAS_rsid, case_alleles, poss_geno, control_alleles, AA, AB, BB):
            genotypes = set()

            all_pairs_set.add((GWAS_rsid, outside_rsid))

            @timeit
            def all_iterations():
                start_time = time()
                for i in range(self.iterations):
                    #logging.debug("this is the {} iteration".format(i))

                    iteration_processor(
                        self, GWAS_rsid, case_alleles, poss_geno, control_alleles, AA, AB, BB, genotypes)

                    #logging.debug("each permutation iteration, time elapsed: {0}".format(exec_time))
                end_time = time()
                exec_time = end_time - start_time
                logging.debug("Performed {} iterations in {} seconds, average is: {} iterations per second".format(
                    self.iterations, exec_time, self.iterations / exec_time))

            def apply_filter():
                all_info_filename = "{}_ALL_INFO_{}_iterations".format(
                    self.p_value_filename, self.iterations)
                # open file here but write in the for loop so you don't
                # open/close files multiple times
                with open(all_info_filename, "a+") as f:
                    labels = ["GWAS_rsID", "outside_rsID", "GWAS_geno", "outside_geno", "case_count", "case_total", "case_odds", "control_count", "control_total", "control_odds", "odds_ratio_GWAS_and_outside", "num_perm_lower",
                              "num_perm_higher", "iterations", "odds_ratio_GWAS", "outside_cause_higher_or_lower_risk", "iter_used_for_p_value", "p_value", "p_value_no_zero", "CI_lower", "CI_higher", "mtc_threshold", "need_more_perm", "status"]
                    f.write("\t".join(labels))
                    f.write("\n")

                    #case_count, case_total, case_odds, control_count, control_total, control_odds, odds_ratio
                    for geno in sorted(genotypes):
                        key = ((GWAS_rsid, outside_rsid), geno)
                        p_value_line = self.p_value_dict[key]
                        # logging.debug("P VALUE LINE IS ", p_value_line)
                        compare_val = p_value_line[-1]
                        odds_ratio = p_value_line[self.odds_ratio_index]

                        if (compare_val != "NA") & (odds_ratio != "NA") & (compare_val != odds_ratio):
                            total = float(p_value_line[self.total_index])
                            # odds_ratio = p_value_line[self.odds_ratio_index]
                            if total != "NA" and total > 0:
                                high_iter = p_value_line[self.high_p_index]
                                low_iter = p_value_line[self.low_p_index]

                                high_p = high_iter / total
                                low_p = low_iter / total

                                real_p_value, iter_used_for_p, outside_cause_higher_or_lower_risk = (
                                    high_p, high_iter, "higher") if odds_ratio > compare_val else (low_p, low_iter, "lower")

                                # the get_compare_info_list method will handle
                                # p_values that are 0s
                                p_value_no_zero, c_i_neg, c_i_pos, mtc_threshold, need_more_perm, status = self.stepwise_filter.get_compare_info_list(
                                    self.iterations, real_p_value, GWAS_rsid, outside_cause_higher_or_lower_risk)

                                write_list = [GWAS_rsid, outside_rsid]
                                write_list.extend(list(geno))
                                write_list += p_value_line
                                write_list += [outside_cause_higher_or_lower_risk, iter_used_for_p, real_p_value, p_value_no_zero, c_i_neg, c_i_pos,
                                               mtc_threshold, need_more_perm, status]
                                write_list = [str(i) for i in write_list]

                                f.write("\t".join(write_list))
                                f.write("\n")

                                if need_more_perm:
                                    need_more_perm_pairs_set.add(
                                        (GWAS_rsid, outside_rsid))
                                    # if we know that this SNP pairs need more permutations, then we don't have to check the other genotypes of this SNP pairs
                                    # if dont need more perm, to find pairs
                                    # that don't need more perm , we still have
                                    # to check for whether the pair is already
                                    # in need_more_perm_pairs_set to prevent
                                    # duplicates, so instead we just do set
                                    # difference: all_pairs -
                                    # need_more_perm_pairs
                                    # break

            all_iterations()
            apply_filter()

            # no_more_perm_set = all
            # self.moreperm_file_maker(need_more_perm_pairs_set, no_more_perm_set)

            if self.output_partial:
                if self.output_odds:
                    self.random_file_maker_partial(
                        outside_rsid, GWAS_rsid, genotypes)

                #self.p_file_maker_partial(outside_rsid, GWAS_rsid, genotypes)
                self.p_file_maker_partial_an(
                    outside_rsid, GWAS_rsid, genotypes, need_more_perm_pairs_set)

        def iteration_processor(self, GWAS_rsid, case_alleles, poss_geno, control_alleles, AA, AB, BB, genotypes):
            # only need case_alleles, control_alleles, GWAS_rsid, outside_rsid,
            # poss_geno
            rand_case, rand_control = randomizer(case_alleles, control_alleles)
            case_count_map = self.single_rsid_combo(
                GWAS_rsid, rand_case, poss_geno, self.case_GWAS)
            control_count_map = self.single_rsid_combo(
                GWAS_rsid, rand_control, poss_geno, self.control_GWAS)

            # def generate_test_data(num_perm, file_name):
            #     logging.debug("in generate_test_data")
            #     labels = ["GWAS_rsid", "outside_rsid", "case_alleles", "control_alleles", "rand_case",
            #               "rand_control", "case_count_map", "control_count_map", "odds_ratio_dict"]
            #     with open(file_name, "a+") as f:
            #         f.write("\t".join(labels))
            #         f.write("\n")
            #         for i in range(num_perm):
            #             logging.debug("generating sample ", i)
            #             rand_case, rand_control = randomizer(
            #                 case_alleles, control_alleles)
            #             case_count_map = self.single_rsid_combo(
            #                 GWAS_rsid, rand_case, poss_geno, self.case_GWAS)
            #             control_count_map = self.single_rsid_combo(
            #                 GWAS_rsid, rand_control, poss_geno, self.control_GWAS)
            #             odds_ratio_dict = self.log_odds_calculator(
            #                 case_count_map, control_count_map)

            #             write_list = [GWAS_rsid, outside_rsid, case_alleles, control_alleles,
            #                           rand_case, rand_control, case_count_map, control_count_map, odds_ratio_dict]
            #             write_list = [str(element) for element in write_list]

            #             f.write("\t".join(write_list))
            #             f.write("\n")

            #     logging.debug("finished writing {} samples".format(num_perm))

            # def generate_test_dicts(folder_name):
            #     odds_ratio_dict = self.log_odds_calculator(
            #         case_count_map, control_count_map)

            #     labels = ["case_alleles", "control_alleles", "rand_case", "rand_control",
            #               "case_count_map", "control_count_map", "odds_ratio_dict", ]

            #     save_dicts_list = [case_alleles, control_alleles, rand_case,
            #                        rand_control, case_count_map, control_count_map, odds_ratio_dict]

            #     # save GO specific dicts to separate folders
            #     self.try_make_folder(folder_name)
            #     for dict_name, dict_to_save in zip(labels, save_dicts_list):
            #         dict_file_name = "{}_{}".format(folder_name, dict_name)
            #         self.save_dict(dict_to_save, folder_name, dict_file_name)
            #     logging.debug("generate_test_dicts at ", folder_name)

            # testing
            #----
            #generate_test_data(10000, "randomized_test_data_GO_{}_{}".format(GWAS_rsid, outside_rsid))
            #generate_test_data(1, "randomized_test_data_all_GO_combo")
            #generate_test_dicts("{}_{}".format(GWAS_rsid, outside_rsid))
            #-----
            if self.odds_ratio_type == 3:
                odds_ratio_dict = self.log_odds_calculator(
                    case_count_map, control_count_map)
                for geno_combo in odds_ratio_dict:
                    odds_ratio = odds_ratio_dict[geno_combo]

                    p_value_key = ((GWAS_rsid, outside_rsid), geno_combo)
                    genotypes.add(geno_combo)

                    if self.output_odds:
                        self.odds_ratio_dict[p_value_key].append(odds_ratio)

                    p_value_line = self.p_value_dict[p_value_key]
                    nonrandom_odds_ratio = p_value_line[self.odds_ratio_index]

                    if odds_ratio != "NA" and nonrandom_odds_ratio != "NA":
                        if odds_ratio <= nonrandom_odds_ratio:
                            p_value_line[self.low_p_index] += 1
                        if odds_ratio >= nonrandom_odds_ratio:
                            p_value_line[self.high_p_index] += 1
                        p_value_line[self.total_index] += 1

            # extract information for linear odds ratio
            elif self.odds_ratio_type == 2:
                del case_count_map["NA_indices"]
                del control_count_map["NA_indices"]

                del case_count_map["total"]
                del control_count_map["total"]

                GWAS_dict_cases = {}
                for key in case_count_map:
                    GWAS, outside = key
                    if GWAS not in GWAS_dict_cases:
                        GWAS_dict_cases[GWAS] = {}
                    GWAS_dict_cases[GWAS][outside] = float(
                        len(case_count_map[key]))

                GWAS_dict_controls = {}
                for key in control_count_map:
                    GWAS, outside = key
                    if GWAS not in GWAS_dict_controls:
                        GWAS_dict_controls[GWAS] = {}
                    GWAS_dict_controls[GWAS][outside] = float(
                        len(control_count_map[key]))

                for GWAS_geno in GWAS_dict_cases:
                    inner_case = GWAS_dict_cases[GWAS_geno]
                    inner_control = GWAS_dict_controls[GWAS_geno]

                    num_A_cases = inner_case[AA] * 2 + inner_case[AB]
                    num_B_cases = inner_case[BB] * 2 + inner_case[AB]

                    num_A_controls = inner_control[AA] * 2 + inner_control[AB]
                    num_B_controls = inner_control[BB] * 2 + inner_control[AB]

                    if num_A_cases + num_A_controls <= num_B_cases + num_B_controls:
                        first_cases, first_controls = num_A_cases, num_A_controls
                        second_cases, second_controls = num_B_cases, num_B_controls
                    else:
                        first_cases, first_controls = num_B_cases, num_B_controls
                        second_cases, second_controls = num_A_cases, num_A_controls

                    odds_ratio = self.odds_ratio_calculator(
                        first_cases, second_cases, first_controls, second_controls)

                    p_value_key = ((GWAS_rsid, outside_rsid), (GWAS_geno,))
                    genotypes.add((GWAS_geno,))

                    if self.output_odds:
                        self.odds_ratio_dict[p_value_key].append(odds_ratio)

                    p_value_line = self.p_value_dict[p_value_key]
                    nonrandom_odds_ratio = p_value_line[self.odds_ratio_index]

                    if odds_ratio != "NA" and nonrandom_odds_ratio != "NA":
                        if odds_ratio <= nonrandom_odds_ratio:
                            p_value_line[self.low_p_index] += 1
                        if odds_ratio >= nonrandom_odds_ratio:
                            p_value_line[self.high_p_index] += 1
                        p_value_line[self.total_index] += 1

            # extract information for nonlinear odds ratio
            else:
                case_total = float(
                    case_count_map["total"] - len(case_count_map["NA_indices"]))
                control_total = float(
                    control_count_map["total"] - len(control_count_map["NA_indices"]))

                del case_count_map["total"]
                del control_count_map["total"]

                del case_count_map["NA_indices"]
                del control_count_map["NA_indices"]

                for geno_combo in case_count_map:
                    case_count = len(case_count_map[geno_combo])
                    control_count = len(control_count_map[geno_combo])

                    odds_ratio = self.odds_ratio_calculator(
                        case_count, control_count, case_total, control_total)

                    p_value_key = ((GWAS_rsid, outside_rsid), geno_combo)
                    genotypes.add(geno_combo)

                    if self.output_odds:
                        self.odds_ratio_dict[
                            p_value_key].append(odds_ratio)

                    p_value_line = self.p_value_dict[p_value_key]
                    nonrandom_odds_ratio = p_value_line[
                        self.odds_ratio_index]

                    if odds_ratio != "NA" and nonrandom_odds_ratio != "NA":
                        if odds_ratio <= nonrandom_odds_ratio:
                            p_value_line[self.low_p_index] += 1
                        if odds_ratio >= nonrandom_odds_ratio:
                            p_value_line[self.high_p_index] += 1
                        p_value_line[self.total_index] += 1
            return genotypes

        def make_report_file():

            report_file_name = "{}_report".format(self.p_value_filename)
            num_filtered_file_name = "{}_num_filtered_per_interval".format(
                self.p_value_filename)

            num_filtered_headers = ["Iterations", "Num no more perm",
                                    "Num more perm", "Percent more perm", "time", "Percent no more perm"]

            if not os.path.isfile(num_filtered_file_name):
                with open(num_filtered_file_name, "w+") as f:
                    f.write("\t".join(num_filtered_headers))
                    f.write("\n")

            # TODO: Make this more legible

            result_info = [self.iterations, len(all_pairs_set), len(no_more_perm_set), len(need_more_perm_pairs_set), (float(len(need_more_perm_pairs_set)) / len(
                self.all_SNP_pairs_in_pairing_file_set)) * 100, time(), (float(len(no_more_perm_set)) / len(self.all_SNP_pairs_in_pairing_file_set)) * 100]

            with open(report_file_name, "a+") as f:
                f.write("Step finished at {5}. Ran {0} iterations, started with {1} pairs, filtered {2} pairs, {3} pairs still need more permutations. Percent filtered since last step: {6} Percent still need more permutations: {4} percent ".format(*result_info
                                                                                                                                                                                                                                                        ))
                f.write("\n")
                f.write("-----------------------")
                f.write("\n")

            with open(num_filtered_file_name, "a+") as f:
                f.write("\t".join(list(map(str, result_info))))
                f.write("\n")

        # these sets get updated in GWAS_rsid_processor function
        need_more_perm_pairs_set = set()
        all_pairs_set = set()

        for outside_rsid in self.case_outside:
            # TODO: make set updating explicit
            #need_more_perm_pairs_set, all_pairs_set = outside_rsid_processor(self, outside_rsid)
            outside_rsid_processor(self, outside_rsid)

        no_more_perm_set = all_pairs_set - need_more_perm_pairs_set
        self.moreperm_file_maker(need_more_perm_pairs_set, no_more_perm_set)

        no_more_perm_outside_rsids = [outside_rsid for (
            _, outside_rsid) in no_more_perm_set]

        # delete outside_rsids that don't need more permutations
        for outside_rsid_no_longer_need in no_more_perm_outside_rsids:
            # , self.control_outside, self.O2G]: ##deleted this because it caused problem in iterations
            for stored_dict in [self.case_outside]:
                stored_dict.pop(outside_rsid_no_longer_need)

        # Disabled this because some pipelines don't have self.all_SNP_pairs_in_pairing_file_set due to different pairing_maker() function.
        # make_report_file()
        return need_more_perm_pairs_set

    def single_rsid_combo(self, GWAS_rsid, rand_case, poss_geno, GWAS_dict):
        """
        performs odds ratio calculation for a single GWAS-outside rsID pair
        """
        geno_combo_to_count = defaultdict(list)
        # logging.debug(GWAS_dict)
        current_GWAS_map = GWAS_dict[GWAS_rsid]
        GWAS_alleles = current_GWAS_map["alleles"]

        geno_combo_to_count["total"] = current_GWAS_map["total"]
        geno_combo_to_count["NA_indices"] = current_GWAS_map["NA_indices"][:]

        for GWAS_allele in GWAS_alleles:
            columns = current_GWAS_map[GWAS_allele]
            for col in columns:
                outside_allele = rand_case[col]
                if outside_allele == self.NA:
                    geno_combo_to_count["NA_indices"].append(col)
                else:
                    new_key = (GWAS_allele, outside_allele)
                    geno_combo_to_count[new_key].append(col)

            for poss_outside in poss_geno:
                check = (GWAS_allele, poss_outside)
                geno_combo_to_count.setdefault(check, [])

        return geno_combo_to_count

    def random_file_maker_partial(self, outside_rsid, GWAS_rsid, genotypes):
        """
        writes partial random odds ratio information to output file
        """
        f = open(self.rand_filename, "a")

        for geno in sorted(genotypes):
            write_list = [GWAS_rsid, outside_rsid]
            write_list.extend(list(geno))

            key = ((GWAS_rsid, outside_rsid), geno)
            for elt in self.odds_ratio_dict[key]:
                if isinstance(elt, float):
                    elt = round(elt, self.num_decimals)
                write_list.append(str(elt))
            f.write("\t".join(write_list))
            f.write("\n")

            del self.odds_ratio_dict[key]

        f.close()

    def random_file_maker(self):
        """
        creates file based off odds_ratio_dict
        """
        f = open(self.rand_filename, "w")

        sorted_key = sorted(self.odds_ratio_dict.keys())
        for key in sorted_key:
            (GWAS_rsid, outside_rsid), geno_key = key
            write_list = [GWAS_rsid, outside_rsid]
            write_list.extend(list(geno_key))

            for elt in self.odds_ratio_dict[key]:
                if isinstance(elt, float):
                    elt = round(elt, self.num_decimals)
                write_list.append(str(elt))
            f.write("\t".join(write_list))
            f.write("\n")

        f.close()


    def p_file_maker_partial_an(self, outside_rsid, GWAS_rsid, genotypes, need_more_perm_pairs_set):
        iteration_filename = "{}_{}_iterations".format(
            self.p_value_filename, str(self.iterations))

        snp_pair = (GWAS_rsid, outside_rsid)
        need_more_perm = snp_pair in need_more_perm_pairs_set

        with open(iteration_filename, "a+") as f:
            for geno in sorted(genotypes):
                write_list = [GWAS_rsid, outside_rsid]
                write_list.extend(list(geno))

                key = ((GWAS_rsid, outside_rsid), geno)
                p_value_line = self.p_value_dict[key]
                total = p_value_line[self.total_index]
                if total != "NA" and total > 0:
                    # calculate fraction
                    low_p = float(p_value_line[self.low_p_index]) / total
                    high_p = float(p_value_line[self.high_p_index]) / total
                else:
                    low_p = "NA"
                    high_p = "NA"

                for i in range(len(p_value_line)):
                    curr = p_value_line[i]
                    if i in self.int_indices and curr != "NA":
                        curr = int(curr)
                    if isinstance(curr, float) and i != self.low_p_index and i != self.high_p_index:
                        curr = round(curr, self.num_decimals)

                    # replace the randomization count by fraction for output
                    # file
                    if i == self.low_p_index:
                        write_list.append(str(low_p))
                    elif i == self.high_p_index:
                        write_list.append(str(high_p))
                    else:
                        write_list.append(str(curr))
                    if i == len(p_value_line) - 1:
                        odds_ratio = p_value_line[self.odds_ratio_index]
                        if curr == "NA" or odds_ratio == "NA":
                            write_list.append("NA")
                        elif curr > odds_ratio:
                            write_list.append((str(low_p)))
                        elif curr < odds_ratio:
                            write_list.append((str(high_p)))
                        else:
                            write_list.append("EQ")

                write_list = [str(i) for i in write_list]
                # GWAS_rsid
                f.write("\t".join(write_list))
                f.write("\n")
                if not need_more_perm:
                    del self.p_value_dict[key]

        logging.debug(f"created file: {iteration_filename}")

    def moreperm_file_maker(self, moreperm_pairs_set, no_more_perm_set):
        moreperm_filename = "{}_moreperm_{}_iterations".format(
            self.p_value_filename, str(self.iterations))
        no_more_perm_filename = "{}_nomoreperm_{}_iterations".format(
            self.p_value_filename, str(self.iterations))

        with open(moreperm_filename, "w+") as f:
            for pair in sorted(moreperm_pairs_set):
                f.write(" ".join([str(ele) for ele in pair]))
                f.write("\n")

        with open(no_more_perm_filename, "w+") as f:
            for pair in sorted(no_more_perm_set):
                f.write(" ".join([str(ele) for ele in pair]))
                f.write("\n")

    def p_file_maker(self):
        """
        creates file based off p_value_dict
        """
        f = open(self.p_value_filename, "w")

        sorted_key = sorted(self.p_value_dict.keys())
        for key in sorted_key:
            (GWAS_rsid, outside_rsid), geno_key = key
            write_list = [GWAS_rsid, outside_rsid]
            write_list.extend(list(geno_key))

            p_value_line = self.p_value_dict[key]
            total = p_value_line[self.total_index]
            if total != "NA" and total > 0:
                p_value_line[self.low_p_index] = float(
                    p_value_line[self.low_p_index]) / total
                p_value_line[self.high_p_index] = float(
                    p_value_line[self.high_p_index]) / total
            for i in range(len(p_value_line)):
                curr = p_value_line[i]
                if i in self.int_indices and curr != "NA":
                    curr = int(curr)
                if isinstance(curr, float) and i != self.low_p_index and i != self.high_p_index:
                    curr = round(curr, self.num_decimals)
                write_list.append(str(curr))
                if i == len(p_value_line) - 1:
                    odds_ratio = p_value_line[self.odds_ratio_index]
                    if curr == "NA" or odds_ratio == "NA":
                        write_list.append("NA")
                    elif curr > odds_ratio:
                        write_list.append(
                            (str(p_value_line[self.low_p_index])))
                    elif curr < odds_ratio:
                        write_list.append(
                            (str(p_value_line[self.high_p_index])))
                    else:
                        write_list.append("EQ")
            f.write("\t".join(write_list))
            f.write("\n")

        f.close()

    @timeit
    def read_input_files(self):

        if self.check_start:
            logging.info("CHECKING FORMAT OF GENOTYPE FILE(S)")
            self.check_format()

        logging.info("READING PAIRING FILE")
        with swap_attr(self, {"delim":"pairing_delim"}):
            self.pairing_maker()
        if not self.one_pair:
            self.original_num_pairs = pd.read_csv(
                self.pair_filename, header=None, sep=" ").shape[0]
        logging.info("READING SAMPLE FILE(S) (if found)")

        with swap_attr(self, {"delim":"sample_delim"}):
            self.read_sample_files()

    # TODO: remove this function from run functions and call it from main
    # program
    @timeit
    def run_pre_random(self):
        start_time = time()
        #case_control_dict_filename = "case_control_dict_{}".format(self.hash)
        case_control_dict_filename = "case_control_dict"
        # try to load file first, else create from scratch and then save for
        # next time

        # TODO: turn this into a decorator: dec(function, *data_needed)  to try
        # loading data_needed from disk from a previous run before function
        # execution
        try:
            # case_control_dict = self.load_dict(case_control_dict_filename, load_dir=self.args_input_folder_path)
            try:
                # try load single pair dict
                case_control_dict = self.load_dict(case_control_dict_filename)
                found_single_pair_dict = True
                logging.info("Loaded single pair dict in rsid folder path")
            except FileNotFoundError:
                logging.info(f"Did not find single pair dict, loading entire dict")
                case_control_dict = self.load_dict(case_control_dict_filename, load_dir=self.args_input_folder_path)
                found_single_pair_dict = False

            self.case_GWAS, self.case_outside = case_control_dict["case"]
            self.control_GWAS, self.control_outside = case_control_dict[
                "control"]

            self.case_GWAS = {key: self.case_GWAS[key] for key in self.GWAS_set}
            self.case_outside = {key: self.case_outside[key] for key in self.outside_set}


            self.control_GWAS, self.control_outside = case_control_dict[
                "control"]

            self.control_GWAS = {key: self.control_GWAS[key] for key in self.GWAS_set}
            self.control_outside = {key: self.control_outside[key] for key in self.outside_set}

            #if did not find single pair dict but found all pairs dict, then just save what's relevant for this pair in the pair directory
            if not found_single_pair_dict:
                logging.info("Outputting single pair dict in rsid folder path")
                case_control_dict = {"case": (self.case_GWAS, self.case_outside),
                                 "control": (self.control_GWAS, self.control_outside)}

                self.save_dict(case_control_dict, self.working_dir,
           case_control_dict_filename)

            logging.info("loaded G_O_dict from file")
        except FileNotFoundError:
            logging.info(f"File for G_O_dict not found in current directory: {os.getcwd()}, creating from scratch")

            with cd(self.input_folder_path):
                logging.info("READING GENOTYPE FILE(S)")
                #swap out delimiter attribute to read genetic file, which might have a different delimiter than sample file or pairing file
                with swap_attr(self, {"delim":"gen_delim"}):

                    self.G_O_dict_maker()


            # # catching concurrent I/O error. If fail, sleep and retry
            # while True:
            #   try:
            #       self.G_O_dict_maker()
            #
            #   except FileNotFoundError :
            #       logging.debug("not found file")
            #       sleep(5)
            #       continue
            #   break

            # saving dict results from G_O_dict_maker

            case_control_dict = {"case": (self.case_GWAS, self.case_outside),
                                 "control": (self.control_GWAS, self.control_outside)}

            self.save_dict(case_control_dict, self.working_dir,
                           case_control_dict_filename)
            # free up space
            del case_control_dict

        self.GWAS_OR_calc()

        # case_combined_dict_filename = "case_combined_dict_{}".format(self.hash)
        # control_combined_dict_filename = "control_combined_dict_{}".format(
        #     self.hash)

        # # TODO: turn this into a decorator: dec(function, *data_needed)  to try
        # # loading data_needed from disk from a previous run before function
        # # execution
        # try:
        #     self.case_combined = self.load_dict(case_combined_dict_filename)
        #     self.control_combined = self.load_dict(
        #         control_combined_dict_filename)
        #     logging.info("Loaded case_combined and control_combined from file")
        # except FileNotFoundError:
        #     logging.info(
        #         "File for case_combined and/or control_combined not found, creating and saving both dicts")
        #     self.combine_dicts()
        #     self.save_dict(self.case_combined, self.working_dir,
        #                    case_combined_dict_filename)
        #     self.save_dict(self.control_combined, self.working_dir,
        #                    control_combined_dict_filename)

        self.combine_dicts()

        end_time = time()
        exec_time = end_time - start_time
        logging.info(
            "work before the randomization, time elapsed: {0}".format(exec_time))


    @timeit
    def run_random(self):
        def output_shared_dicts():
            shared_dicts = [self.case_GWAS, self.control_GWAS,
                            self.NA, self.covariate_data]
            shared_dicts_labels = [
                "self_case_GWAS", "self_control_GWAS", "self_NA", "self_covariate_data"]

            # save shared_dicts to outer folder
            for (dict_name, dict_to_save) in zip(shared_dicts_labels, shared_dicts):
                self.save_dict(dict_to_save, self.working_dir, dict_name)

        logging.info("OUTPUTTING SHARED DICTS")
        # testing to pull out just the info we need
        # output_shared_dicts()
        logging.info("PERFORMING RANDOMIZATIONS")
#       self.perform_randomizations()

        # testing
        # move to working directory:
        # for i in [1000, 10_000, 100_000, 1_000_000, 10_000_000, 100_000_000]:  # + ([100] * 1000):
        #   self.iterations = i
        more_perm = self.perform_randomizations_an()
        if len(more_perm) == 0:
            logging.debug("no more snp that needs permutation beyond {} iterations".format(
                self.iterations))

        start_time = time()
        logging.info("CREATING OUTPUT FILES")
        if not self.output_partial:
            if self.output_odds:
                self.random_file_maker()
            self.p_file_maker()

        end_time = time()
        exec_time = end_time - start_time
        logging.info(
            "work after the randomization, time elapsed: {0}".format(exec_time))
        logging.info("\nDONE\n")

    def run_with_permutation(self):
        def run_func():
            self.run_pre_random()
            logging.info("Calculating odds ratio for the unpermutated dataset")
            self.nonrand_dict_maker()
            self.run_random()

            with atomic_write("finished") as f:
                f.write("marker to tell orchestrator that this job has finished")

        return timeit(run_func, name="run_func.run_with_permutation")

        # This run version should ONLY be called with ALL SNP pairs for each GWAS SNP to calculate accurate MTC
    def run_create_mtc_table(self, num_case_impute_filtering, num_control_impute_filtering, case_gen, control_gen):
        def run_func():
            # makes self.case_combined and self.control_combined
            make_mtc_output_folder = "input_folder_get_mtc"
            os.makedirs(make_mtc_output_folder)
            with cd(make_mtc_output_folder):
                self.run_pre_random()

                self.filter_df.to_csv(f"filter_df_found_not_found_{self.p_value_filename}", sep="\t",index=False)

                all_snp_pairs_mtc_df, found_pairs, after_filter_1, after_filter_2, filtered_after_1, filtered_after_2 = MtcTable.make_mtc_table_from_dict(
                    self.case_combined, self.control_combined, num_case_impute_filtering, num_control_impute_filtering)
                all_snp_pairs_mtc_df.to_csv("MTC_table_{}".format(
                    self.p_value_filename), index=False, sep='\t')
                after_filter_1 = after_filter_1.sort_values(
                    by=["GWAS_rsid", "outside_rsid"])
                after_filter_2 = after_filter_2.sort_values(
                    by=["GWAS_rsid", "outside_rsid"])

                logging.info("CREATING MTC TABLE")
                with atomic_write("make_mtc_report") as f:
                    f.write("Original number of SNP pairs: {} \n".format(
                        self.original_num_pairs))
                    f.write(
                        "Number of pairs found in genetic files: {} \n".format(found_pairs))
                    f.write("Number of pairs after Filter 1: {} \n".format(
                        after_filter_1[["GWAS_rsid", "outside_rsid"]].drop_duplicates().shape[0]))
                    f.write("Number of pairs after Filter 2: {} \n".format(
                        after_filter_2[["GWAS_rsid", "outside_rsid"]].drop_duplicates().shape[0]))
                after_filter_1.to_csv("After_filter_1_{}".format(
                    self.p_value_filename), index=False, sep=' ')
                after_filter_2.to_csv("After_filter_2_{}".format(
                    self.p_value_filename), index=False, sep=' ')

                filtered_snp_pairs = after_filter_2[
                    ["GWAS_rsid", "outside_rsid"]].drop_duplicates()

                filtered_snp_pairs.to_csv("Filtered_SNP_pairs_{}_file".format(
                    self.p_value_filename), header=False, index=False, sep=self.pairing_delim)

                filtered_after_1.to_csv(f"Filtered_OUT_AFTER_FILTER1_{self.p_value_filename}_file",index=False, sep="\t")
                filtered_after_2.to_csv(f"Filtered_OUT_AFTER_FILTER2_{self.p_value_filename}_file",index=False, sep="\t")

            logging.info(
                "CREATED MTC TABLE. NOW YOU CAN RUN THE SCRIPT WITH INIT FILE THAT GETS MTC TABLE YOU JUST CREATED")
        return timeit(run_func, name="run_func.run_create_mtc_table")

    def __str__(self):
        rep_dict = {
            "CASE": self.case_filename,
            "CONTROL": self.control_filename,
            "PAIR": self.pair_filename,
            "RAND_OUTPUT": self.rand_filename,
            "P_VAL_OUTPUT": self.p_value_filename,
            "GENO_START": self.col_cut,
            "ITERATIONS": self.iterations,
            "TRIPLES": self.is_triples,
            "RSID_COL": self.rsid_col,
            "SNP_COL": self.snp_col,
            "DELIM": self.delim,
            "OUTPUT_ODDS": self.output_odds
        }
        return str(rep_dict)

    def __repr__(self):
        return self.__str__()



class Case_Control_Pipeline(Pipeline):

    def __init__(self, case, control, pair, rand, p, gs, itr, typ, rsid, snp, delim, or_calc, skip, path, partial, check, case_sample, control_sample, covs, high, low, NA):
        Pipeline.__init__(self, pair, rand, p, gs, itr, typ, rsid, snp,
                          delim, or_calc, skip, path, partial, check, covs, high, low, NA)
        self.case_filename = path + case
        self.control_filename = path + control
        self.case_sample_file = case_sample
        self.control_sample_file = control_sample

    def check_format(self):
        """
        checks format of case and control files
        """
        self.check_format_single_file(self.case_filename)
        self.check_format_single_file(self.control_filename)

    @timeit
    def read_sample_files(self):
        if self.case_sample_file != None and self.control_sample_file != None:
            case_cov_data = self.read_single_sample_file(self.case_sample_file)
            control_cov_data = self.read_single_sample_file(
                self.control_sample_file)
            self.covariate_data = np.vstack((case_cov_data, control_cov_data))

    def read_single_sample_file(self, sample_file):
        f = open(sample_file, "r")

        header = f.readline()
        cov_indices = self.extract_cov_info(header, sample_file)

        f.readline()

        sample_covs = []
        for line in f:
            new_cov_line = []
            for index in cov_indices:
                covariate_data = float(self.column_getter(line, index + 1))
                new_cov_line.append(covariate_data)
            sample_covs.append(new_cov_line)
        sample_covs = np.array(sample_covs)

        f.close()

        return sample_covs

    @timeit
    def G_O_dict_maker(self):
        """
        reads case and control files and creates genetic info data structures
        """
        self.case_GWAS, self.case_outside = self.G_O_dict_maker_single(
            self.case_filename)
        self.control_GWAS, self.control_outside = self.G_O_dict_maker_single(
            self.control_filename)

        if not all([self.case_GWAS, self.case_outside, self.control_GWAS, self.control_outside]):
            raise LookupError(
                "Check messages above. None of the GWAS/outside rsIDs were found in the case/control genotype file.")

        return self.case_GWAS, self.case_outside, self.control_GWAS, self.control_outside

    def G_O_dict_maker_single(self, filename):
        """
        creates GWAS and outside genetic info data structures for a single file
        """
        file = open(filename, "r")

        found_GWAS = set()
        found_outside = set()
        GWAS_map = {}
        outside_map = {}

        line_num = 1

        for i,line in tqdm(enumerate(file)):
            rsid = self.column_getter(line, self.rsid_col)
            if rsid in self.GWAS_set or rsid in self.outside_set:
                pos_1 = self.column_getter(line, self.snp_col)
                if "/" in pos_1:
                    pos_2 = pos_1[2]
                    pos_1 = pos_1[0]
                else:
                    pos_2 = self.column_getter(line, self.snp_col + 1)
                if len(pos_1) > 1 or len(pos_2) > 1:
                    continue

                poss_geno = possible_genotypes(pos_1, pos_2)
                head, alleles = self.column_cutter(line, self.col_cut)
                allele_list = alleles.split(self.delim)
                if self.is_triples:
                    if len(allele_list) % 3 != 0:
                        raise ValueError(
                            "Missing data. Incomplete genotype triple found in line " + str(line_num) + ".")
                    else:
                        allele_list = self.triplicate_converter(
                            allele_list, pos_1, pos_2)

                allele_list = list(map(self.clean_item_cached, allele_list))
                # remove list here
                total = len(list(allele_list))

                if rsid in self.GWAS_set:
                    found_GWAS.add(rsid)
                    new_map = defaultdict(list)
                    new_map["alleles"] = poss_geno
                    new_map["total"] = total
                    for i in range(len(list(allele_list))):
                        curr_allele = allele_list[i]
                        if curr_allele != self.NA:
                            new_map[curr_allele].append(i)
                        else:
                            new_map["NA_indices"].append(i)
                    GWAS_map[rsid] = new_map
                else:
                    found_outside.add(rsid)
                    outside_map[rsid] = allele_list, poss_geno
            line_num += 1
        file.close()

        logging.info(str(len(found_GWAS)) + " out of " + str(len(self.GWAS_set)
                                                             ) + " GWAS rsIDs found in " + filename + ".")
        not_found_GWAS = self.GWAS_set.difference(found_GWAS)
        if len(not_found_GWAS) > 0:
            logging.info("WARNING: The following GWAS rsIDs are found in the pairing file, but not " +
                         filename + ".\n" + str(not_found_GWAS))
            if not_found_GWAS == self.GWAS_set:
                return None, None

        logging.info(str(len(found_outside)) + " out of " + str(len(self.outside_set)
                                                                ) + " outside rsIDs found in " + filename + ".")
        not_found_outside = self.outside_set.difference(found_outside)
        if len(not_found_outside):
            logging.info("WARNING: The following outside rsIDs are found in the pairing file, but not " +
                         filename + ".\n" + str(not_found_outside))
            if not_found_outside == self.outside_set:
                return None, None

        df_data = [*zip(not_found_GWAS, ["not_found"] * len(not_found_GWAS), ["GWAS"]* len(not_found_GWAS)),
            *zip(not_found_outside, ["not_found"] * len(not_found_outside), ["outside"]* len(not_found_outside)),]

        self.filter_df = pd.DataFrame(df_data, columns=["rsID", "reason_for_filtering", "SNP_type"])

        return GWAS_map, outside_map


class Combined_Pipeline(Pipeline):

    def __init__(self, geno, sample, pair, rand, p, gs, itr, typ, rsid, snp, delim, or_calc, skip, path, partial, check, covs, high, low, NA):
        Pipeline.__init__(self, pair, rand, p, gs, itr, typ, rsid, snp,
                          delim, or_calc, skip, path, partial, check, covs, high, low, NA)
        self.geno_filename = path + geno
        self.sample_filename = sample

        self.case_keyword = "case"
        self.exclude_keyword = "missing"

    def check_format(self):
        """
        checks format of .gen file
        """
        self.check_format_single_file(self.geno_filename)

    @timeit
    def read_sample_files(self):
        """
        reads sample file
        returns dictionary mapping line number to case/control int
        returns int representing control status
        """
        out_dict = {}
        f = open(self.sample_filename, "r")

        # assumes first line contains header information
        first_line = f.readline()
        first_line = first_line.strip()
        index_to_find = self.index_getter(first_line, self.case_keyword)
        exclude_index = self.index_getter(first_line, self.exclude_keyword)
        check_exclude = True if exclude_index > -1 else False

        cov_indices = self.extract_cov_info(first_line, self.sample_filename)

        # assumes second line does not contain pertinent information
        trash = f.readline()

        num_set = set()
        cc_count = defaultdict(int)
        line_num = 0
        num_excluded = 0

        sample_covs = []
        for line in f:
            include = True
            if check_exclude:
                include = not int(self.column_getter(line, exclude_index + 1))
            if include:
                case_control = int(self.column_getter(line, index_to_find + 1))
                num_set.add(case_control)
                cc_count[case_control] += 1
                out_dict[line_num] = case_control

                new_cov_line = []
                for index in cov_indices:
                    covariate_data = float(self.column_getter(line, index + 1))
                    new_cov_line.append(covariate_data)
                sample_covs.append(new_cov_line)
            else:
                num_excluded += 1
            line_num += 1
        self.covariate_data = np.array(sample_covs)

        f.close()

        if num_excluded == 0:
            logging.info("ALL SAMPLES INCLUDED")
        else:
            logging.info(str(num_excluded) + " SAMPLES EXCLUDED")

        control_int = min(num_set)

        num_control = cc_count[control_int]
        num_cases = cc_count[control_int + 1]
        logging.info(str(num_control) + " CONTROLS INCLUDED")
        logging.info(str(num_cases) + " CASES INCLUDED")

        self.cc_dict = out_dict
        self.control_int = control_int
        return self.cc_dict, self.control_int
        #return num_control, num_cases

    def cc_splitter(self, allele_list):
        """
        takes in case_control_dict, allele_list (in letter format), and control int
        sorts the allele_list into control and case groups based index's value in case_control_dict
        returns two lists: list of control alleles and list of case alleles
        """
        swap_dict = {self.control_int: 0, self.control_int + 1: 1}
        control_list = []
        case_list = []
        list_of_lists = [control_list, case_list]
        for i in range(len(allele_list)):
            if i in self.cc_dict:
                inner_ind = self.cc_dict[i]
                ind = swap_dict[inner_ind]
                list_of_lists[ind].append(allele_list[i])
        return control_list, case_list

    @timeit
    def G_O_dict_maker(self):
        """
        reads .gen file and creates genetic info data structures
        """
        self.case_GWAS, self.case_outside, self.control_GWAS, self.control_outside = self.G_O_dict_maker_single()
        if not all([self.case_GWAS, self.case_outside, self.control_GWAS, self.control_outside]):
            raise LookupError(
                "Check messages above. None of the GWAS/outside rsIDs were found in the case/control genotype file.")
        return self.case_GWAS, self.case_outside, self.control_GWAS, self.control_outside

    def G_O_dict_maker_single(self):
        """
        creates case and control GWAS and rsID genetic info data structures
        """
        file = open(self.geno_filename, "r")

        found_GWAS = set()
        found_outside = set()
        found_but_wrong_format = set()

        case_GWAS_map = {}
        case_outside_map = {}

        control_GWAS_map = {}
        control_outside_map = {}

        line_num = 1
        for i,line in tqdm(enumerate(file)):
            rsid = self.column_getter(line, self.rsid_col)
            if rsid in self.GWAS_set or rsid in self.outside_set:
                pos_1 = self.column_getter(line, self.snp_col)
                if "/" in pos_1:
                    pos_2 = pos_1[2]
                    pos_1 = pos_1[0]
                else:
                    pos_2 = self.column_getter(line, self.snp_col + 1)
                if len(pos_1) > 1 or len(pos_2) > 1:
                    found_but_wrong_format.add(rsid)
                    continue

                poss_geno = possible_genotypes(pos_1, pos_2)
                head, alleles = self.column_cutter(line, self.col_cut)
                allele_list = alleles.split(self.delim)
                if self.is_triples:
                    if len(allele_list) % 3 != 0:
                        raise ValueError(
                            "Missing data. Incomplete genotype triple found in line " + str(line_num) + ".")
                    else:
                        allele_list = self.triplicate_converter(
                            allele_list, pos_1, pos_2)

                allele_list = list(map(self.clean_item_cached, allele_list))
                control_allele_list, case_allele_list = self.cc_splitter(
                    allele_list)

                case_total = len(case_allele_list)
                control_total = len(control_allele_list)

                if rsid in self.GWAS_set:
                    found_GWAS.add(rsid)
                    new_case_map = defaultdict(list)
                    new_case_map["alleles"] = poss_geno
                    new_case_map["total"] = case_total
                    for i in range(len(case_allele_list)):
                        curr_allele = case_allele_list[i]
                        if curr_allele != self.NA:
                            new_case_map[curr_allele].append(i)
                        else:
                            new_case_map["NA_indices"].append(i)

                    new_control_map = defaultdict(list)
                    new_control_map["alleles"] = poss_geno
                    new_control_map["total"] = control_total
                    for i in range(len(control_allele_list)):
                        curr_allele = control_allele_list[i]
                        if curr_allele != self.NA:
                            new_control_map[curr_allele].append(i)
                        else:
                            new_control_map["NA_indices"].append(i)

                    case_GWAS_map[rsid] = new_case_map
                    control_GWAS_map[rsid] = new_control_map
                else:
                    found_outside.add(rsid)
                    case_outside_map[rsid] = case_allele_list, poss_geno
                    control_outside_map[rsid] = control_allele_list, poss_geno
            line_num += 1
            if found_GWAS == self.GWAS_set and found_outside == self.outside_set:
                break
        file.close()

        logging.info(str(len(found_GWAS)) + " out of " + str(len(self.GWAS_set)
                                                             ) + " GWAS rsIDs found in " + self.geno_filename + ".")
        found_but_wrong_format_GWAS = found_but_wrong_format & self.GWAS_set

        not_found_GWAS = (self.GWAS_set.difference(found_GWAS)) - found_but_wrong_format_GWAS

        if len(not_found_GWAS) > 0:
            logging.info("WARNING: The following GWAS rsIDs are found in the pairing file, but not " +
                         self.geno_filename + ".\n" + str(not_found_GWAS))

        if len(found_but_wrong_format_GWAS) > 0:
            logging.info(f"WARNING: The following GWAS rsIDs are found but has more than two alleles {found_but_wrong_format_GWAS}")

        if (not_found_GWAS | found_but_wrong_format_GWAS)== self.GWAS_set:
            return None, None, None, None

        logging.info(str(len(found_outside)) + " out of " + str(len(self.outside_set)
                                                                ) + " outside rsIDs found in " + self.geno_filename + ".")
        found_but_wrong_format_outside = found_but_wrong_format & self.outside_set

        not_found_outside = self.outside_set.difference(found_outside) - found_but_wrong_format_outside

        if len(not_found_outside) > 0:
            logging.info("WARNING: The following outside rsIDs are found in the pairing file, but not " +
                         self.geno_filename + ".\n" + str(not_found_outside))

        if found_but_wrong_format_outside:
            logging.info(f"WARNING: The following outside rsIDs are found but has more than 2 alleles {found_but_wrong_format_outside}")


        if (not_found_outside | found_but_wrong_format_outside) == self.outside_set:
                return None, None, None, None

        df_data = [*zip(not_found_GWAS, ["not_found"] * len(not_found_GWAS), ["GWAS"]* len(not_found_GWAS)),
        *zip(found_but_wrong_format_GWAS, ["found_but_wrong_format"] * len(found_but_wrong_format_GWAS), ["GWAS"]* len(found_but_wrong_format_GWAS)),
        *zip(not_found_outside, ["not_found"] * len(not_found_outside), ["outside"]* len(not_found_outside)),
        *zip(found_but_wrong_format_outside, ["found_but_wrong_format"] * len(found_but_wrong_format_outside), ["outside"]* len(found_but_wrong_format_outside))]
        #breakpoint()
        self.filter_df = pd.DataFrame(df_data, columns=["rsID", "reason_for_filtering", "SNP_type"])

        return case_GWAS_map, case_outside_map, control_GWAS_map, control_outside_map

    def __str__(self):
        rep_dict = {
            "GEN": self.geno_filename,
            "SAMPLE": self.sample_filename,
            "PAIR": self.pair_filename,
            "RAND_OUTPUT": self.rand_filename,
            "P_VAL_OUTPUT": self.p_value_filename,
            "GENO_START": self.col_cut,
            "ITERATIONS": self.iterations,
            "TRIPLES": self.is_triples,
            "RSID_COL": self.rsid_col,
            "SNP_COL": self.snp_col,
            "DELIM": self.delim,
            "OUTPUT_ODDS": self.output_odds
        }
        return str(rep_dict)


class Trans_Pipeline(Pipeline):
    """
    ABSTRACT CLASS: use Trans_Combined_Pipeline or Trans_Case_Control_Pipeline instead
    """

    def check_format(self):
        """
        checks format of each .gen file
        """
        filenames = self.get_geno_filenames()
        for file in filenames:
            self.check_format_single_file(file)

    def get_geno_filenames(self):
        """
        extracts filenames from filename files
        """
        geno_filenames = []
        for file_file in self.all_files:
            f = open(file_file, "r")
            new_filenames = [self.filepath + line.strip() for line in f]
            geno_filenames.extend(new_filenames)
            f.close()
        return geno_filenames

    def random_file_maker_partial(self, outside_loc, GWAS_loc, genotypes):
        """
        writes partial random odds ratio information to output file
        """
        f = open(self.rand_filename, "a")

        chr2, outside_rsid = outside_loc
        chr1, GWAS_rsid = GWAS_loc
        for geno in sorted(genotypes):
            write_list = [chr1, GWAS_rsid, chr2, outside_rsid]
            write_list.extend(list(geno))

            key = ((GWAS_loc, outside_loc), geno)
            for elt in self.odds_ratio_dict[key]:
                if isinstance(elt, float):
                    elt = round(elt, self.num_decimals)
                write_list.append(str(elt))
            f.write("\t".join(write_list))
            f.write("\n")

            del self.odds_ratio_dict[key]

        f.close()

    def random_file_maker(self):
        """
        creates file based off odds_ratio_dict
        """
        f = open(self.rand_filename, "w")

        sorted_key = sorted(self.odds_ratio_dict.keys())
        for key in sorted_key:
            ((chr1, GWAS_rsid), (chr2, outside_rsid)), geno_key = key
            write_list = [chr1, GWAS_rsid, chr2, outside_rsid]
            write_list.extend(list(geno_key))
            for elt in self.odds_ratio_dict[key]:
                if isinstance(elt, float):
                    elt = round(elt, self.num_decimals)
                write_list.append(str(elt))
            f.write("\t".join(write_list))
            f.write("\n")

        f.close()

    def p_file_maker_partial(self, outside_loc, GWAS_loc, genotypes):
        """
        writes partial p_value information to output file
        """
        # TODO: turn this into a property shared by all pipelines
        f = open("{}_{}_iterations".format(
            self.p_value_filename, str(self.iterations)), "a")

        chr2, outside_rsid = outside_loc
        chr1, GWAS_rsid = GWAS_loc
        for geno in genotypes:
            write_list = [chr1, GWAS_rsid, chr2, outside_rsid]
            write_list.extend(list(geno))

            key = ((GWAS_loc, outside_loc), geno)
            p_value_line = self.p_value_dict[key]
            total = p_value_line[self.total_index]
            if total != "NA" and total > 0:
                p_value_line[self.low_p_index] = float(
                    p_value_line[self.low_p_index]) / total
                p_value_line[self.high_p_index] = float(
                    p_value_line[self.high_p_index]) / total
            for i in range(len(p_value_line)):
                curr = p_value_line[i]
                if i in self.int_indices and curr != "NA":
                    curr = int(curr)
                if isinstance(curr, float) and i != self.low_p_index and i != self.high_p_index:
                    curr = round(curr, self.num_decimals)
                write_list.append(str(curr))
                if i == len(p_value_line) - 1:
                    odds_ratio = p_value_line[self.odds_ratio_index]
                    if curr == "NA" or odds_ratio == "NA":
                        write_list.append("NA")
                    elif curr > odds_ratio:
                        write_list.append(
                            (str(p_value_line[self.low_p_index])))
                    elif curr < odds_ratio:
                        write_list.append(
                            (str(p_value_line[self.high_p_index])))
                    else:
                        write_list.append("EQ")
            f.write("\t".join(write_list))
            f.write("\n")

            del self.p_value_dict[key]

        f.close()

    def p_file_maker(self):
        """
        creates file based off p_value_dict
        """
        f = open(self.p_value_filename, "w")

        sorted_key = sorted(self.p_value_dict.keys())
        for key in sorted_key:
            ((chr1, GWAS_rsid), (chr2, outside_rsid)), geno_key = key
            write_list = [chr1, GWAS_rsid, chr2, outside_rsid]
            write_list.extend(list(geno_key))

            p_value_line = self.p_value_dict[key]
            total = p_value_line[self.total_index]
            if total != "NA" and total > 0:
                p_value_line[self.low_p_index] = float(
                    p_value_line[self.low_p_index]) / total
                p_value_line[self.high_p_index] = float(
                    p_value_line[self.high_p_index]) / total
            for i in range(len(p_value_line)):
                curr = p_value_line[i]
                if i in self.int_indices and curr != "NA":
                    curr = int(curr)
                if isinstance(curr, float) and i != self.low_p_index and i != self.high_p_index:
                    curr = round(curr, self.num_decimals)
                write_list.append(str(curr))
                if i == len(p_value_line) - 1:
                    odds_ratio = p_value_line[self.odds_ratio_index]
                    if curr == "NA" or odds_ratio == "NA":
                        write_list.append("NA")
                    elif curr > odds_ratio:
                        write_list.append(
                            (str(p_value_line[self.low_p_index])))
                    elif curr < odds_ratio:
                        write_list.append(
                            (str(p_value_line[self.high_p_index])))
                    else:
                        write_list.append("EQ")
            f.write("\t".join(write_list))
            f.write("\n")

        f.close()


class Trans_Combined_Pipeline(Combined_Pipeline, Trans_Pipeline):

    def __init__(self, file_file, sample, pair, rand, p, gs, itr, typ, rsid, snp, delim, or_calc, skip, path, partial, check, covs, high, low, NA):
        Combined_Pipeline.__init__(self, file_file, sample, pair, rand, p, gs, itr, typ,
                                   rsid, snp, delim, or_calc, skip, path, partial, check, covs, high, low, NA)
        self.file_file = file_file
        self.all_files = [file_file]

    def check_format(self):
        Trans_Pipeline.check_format()

    @timeit
    def pairing_maker(self):
        """
        reads pairing file and creates GWAS=>outside rsID dictionary,
        outside=>GWAS rsID dictionary, GWAS rsID set, outside rsID set,
        and rsID=>chr dictionary
        returns these five data structures
        """

        self.G2O = defaultdict(list)
        self.O2G = defaultdict(list)
        self.GWAS_set = set()
        self.outside_set = set()

        self.file_set = set()

        with open_file_or_string(self.pair_filename, "r") as pairing_file:
            for line in pairing_file:
                rsA, geno1, rsB, geno2 = line.strip().split(self.delim)

                absolute_GWAS_path = self.filepath + geno1
                absolute_outside_path = self.filepath + geno2

                GWAS_key = (absolute_GWAS_path, rsA)
                outside_key = (absolute_outside_path, rsB)

                self.GWAS_set.add(GWAS_key)
                self.outside_set.add(outside_key)

                self.G2O[GWAS_key].append(outside_key)
                self.O2G[outside_key].append(GWAS_key)

                self.file_set.add(absolute_GWAS_path)
                self.file_set.add(absolute_outside_path)

        if len(self.G2O) == 0 or len(self.O2G) == 0:
            raise IOError("WARNING: SNP pair file is EMPTY")

        return self.G2O, self.O2G, self.GWAS_set, self.outside_set

    @timeit
    def G_O_dict_maker(self):
        """
        creates case and control GWAS and outside genetic info data structures
        reads through all .gen files, updating data structures after every pass
        """
        self.case_GWAS = {}
        self.case_outside = {}
        self.control_GWAS = {}
        self.control_outside = {}

        found_GWAS = set()
        found_outside = set()

        gen_filenames = self.get_geno_filenames()

        skipped_files = []

        for file in gen_filenames:
            if file in self.file_set:
                new_case_GWAS, new_case_outside, new_control_GWAS, new_control_outside, new_found_GWAS, new_found_outside = self.G_O_dict_maker_single(
                    file)
                self.case_GWAS.update(new_case_GWAS)
                self.case_outside.update(new_case_outside)
                self.control_GWAS.update(new_control_GWAS)
                self.control_outside.update(new_control_outside)
                found_GWAS.update(new_found_GWAS)
                found_outside.update(new_found_outside)
            else:
                skipped_files.append(file)

        if len(skipped_files) > 0:
            logging.info(
                "WARNING: The following files were skipped because they do not appear in the pairing file.\n" + str(skipped_files))

        error = False
        logging.info(str(len(found_GWAS)) + " out of " +
                     str(len(self.GWAS_set)) + " GWAS rsIDs found.")
        not_found_GWAS = self.GWAS_set.difference(found_GWAS)
        if len(not_found_GWAS) > 0:
            logging.info(
                "WARNING: The following GWAS rsIDs are found in the pairing file, but not in their genotype file.\n" + str(not_found_GWAS))
            if not_found_GWAS == self.GWAS_set:
                error = True

        logging.info(str(len(found_outside)) + " out of " +
                     str(len(self.outside_set)) + " outside rsIDs found.")
        not_found_outside = self.outside_set.difference(found_outside)
        if len(not_found_outside):
            logging.info(
                "WARNING: The following outside rsIDs are found in the pairing file, but in their genotype file.\n" + str(not_found_outside))
            if not_found_outside == self.outside_set:
                error = True

        if error:
            raise LookupError(
                "Check messages above. None of the GWAS/outside rsIDs were found in the .gen files.")

        return self.case_GWAS, self.case_outside, self.control_GWAS, self.control_outside

    def G_O_dict_maker_single(self, filename):
        """
        creates intermediate data structures for every run through a .gen file
        """
        file = open(filename, "r")

        found_GWAS = set()
        found_outside = set()

        case_GWAS_map = {}
        case_outside_map = {}

        control_GWAS_map = {}
        control_outside_map = {}

        line_num = 1
        for line in file:
            rsid = self.column_getter(line, self.rsid_col)
            check_key = (filename, rsid)
            if check_key in self.GWAS_set or check_key in self.outside_set:
                pos_1 = self.column_getter(line, self.snp_col)
                if "/" in pos_1:
                    pos_2 = pos_1[2]
                    pos_1 = pos_1[0]
                else:
                    pos_2 = self.column_getter(line, self.snp_col + 1)
                if len(pos_1) > 1 or len(pos_2) > 1:
                    continue

                poss_geno = possible_genotypes(pos_1, pos_2)
                head, alleles = self.column_cutter(line, self.col_cut)
                allele_list = alleles.split(self.delim)
                if self.is_triples:
                    if len(allele_list) % 3 != 0:
                        raise ValueError(
                            "Missing data. Incomplete genotype triple found in line " + str(line_num) + ".")
                    else:
                        allele_list = self.triplicate_converter(
                            allele_list, pos_1, pos_2)

                allele_list = list(map(self.clean_item_cached, allele_list))
                control_allele_list, case_allele_list = self.cc_splitter(
                    allele_list)

                case_total = len(case_allele_list)
                control_total = len(control_allele_list)

                if check_key in self.GWAS_set:
                    found_GWAS.add(check_key)
                    new_case_map = defaultdict(list)
                    new_case_map["alleles"] = poss_geno
                    new_case_map["total"] = case_total
                    for i in range(len(case_allele_list)):
                        curr_allele = case_allele_list[i]
                        if curr_allele != self.NA:
                            new_case_map[curr_allele].append(i)
                        else:
                            new_case_map["NA_indices"].append(i)

                    new_control_map = defaultdict(list)
                    new_control_map["alleles"] = poss_geno
                    new_control_map["total"] = control_total
                    for i in range(len(control_allele_list)):
                        curr_allele = control_allele_list[i]
                        if curr_allele != self.NA:
                            new_control_map[curr_allele].append(i)
                        else:
                            new_control_map["NA_indices"].append(i)

                    case_GWAS_map[check_key] = new_case_map
                    control_GWAS_map[check_key] = new_control_map
                else:
                    found_outside.add(check_key)
                    case_outside_map[check_key] = case_allele_list, poss_geno
                    control_outside_map[
                        check_key] = control_allele_list, poss_geno
            line_num += 1
        file.close()

        return case_GWAS_map, case_outside_map, control_GWAS_map, control_outside_map, found_GWAS, found_outside

    def random_file_maker(self):
        Trans_Pipeline.random_file_maker(self)

    def random_file_maker_partial(self, outside_loc, GWAS_loc, genotypes):
        Trans_Pipeline.random_file_maker_partial(
            self, outside_loc, GWAS_loc, genotypes)

    def p_file_maker(self):
        Trans_Pipeline.p_file_maker(self)

    def p_file_maker_partial(self, outside_loc, GWAS_loc, genotypes):
        Trans_Pipeline.p_file_maker_partial(
            self, outside_loc, GWAS_loc, genotypes)


class Trans_Case_Control_Pipeline(Case_Control_Pipeline, Trans_Pipeline):

    def __init__(self, case_file_file, control_file_file, pair, rand, p, gs, itr, typ, rsid, snp, delim, or_calc, skip, path, partial, check, case_sample, control_sample, covs, high, low, NA):
        Case_Control_Pipeline.__init__(self, case_file_file, control_file_file, pair, rand, p, gs, itr, typ, rsid,
                                       snp, delim, or_calc, skip, path, partial, check, case_sample, control_sample, covs, high, low, NA)
        self.case_file_file = case_file_file
        self.control_file_file = control_file_file
        self.all_files = [case_file_file, control_file_file]
        self.case_sample_file = case_sample
        self.control_sample_file = control_sample

    def check_format(self):
        Trans_Pipeline.check_format()

    @timeit
    def pairing_maker(self):
        """
        reads pairing file and creates GWAS=>outside rsID dictionary,
        outside=>GWAS rsID dictionary, GWAS rsID set, outside rsID set,
        and rsID=>chr dictionary
        returns these five data structures
        """
        self.G2O = defaultdict(list)
        self.O2G = defaultdict(list)

        self.case_GWAS_set = set()
        self.case_outside_set = set()

        self.control_GWAS_set = set()
        self.control_outside_set = set()

        self.file_set = set()

        self.swap_key_dict = {}

        with open_file_or_string(self.pair_filename, "r") as pairing_file:
            for line in pairing_file:
                rsA, case1, control1, rsB, case2, control2 = line.strip().split(self.delim)

                abs_GWAS_case_path = self.filepath + case1
                abs_GWAS_control_path = self.filepath + control1

                abs_outside_case_path = self.filepath + case2
                abs_outside_control_path = self.filepath + control2

                GWAS_case_key = (abs_GWAS_case_path, rsA)
                GWAS_control_key = (abs_GWAS_control_path, rsA)

                outside_case_key = (abs_outside_case_path, rsB)
                outside_control_key = (abs_outside_control_path, rsB)

                self.case_GWAS_set.add(GWAS_case_key)
                self.case_outside_set.add(outside_case_key)

                self.control_GWAS_set.add(GWAS_control_key)
                self.control_outside_set.add(outside_control_key)

                self.file_set.add(abs_GWAS_case_path)
                self.file_set.add(abs_GWAS_control_path)
                self.file_set.add(abs_outside_case_path)
                self.file_set.add(abs_outside_control_path)

                GWAS_key = (case1 + "," + control1, rsA)
                outside_key = (case2 + "," + control2, rsB)

                self.swap_key_dict[GWAS_case_key] = GWAS_key
                self.swap_key_dict[GWAS_control_key] = GWAS_key

                self.swap_key_dict[outside_case_key] = outside_key
                self.swap_key_dict[outside_control_key] = outside_key

                self.G2O[GWAS_key].append(outside_key)
                self.O2G[outside_key].append(GWAS_key)

        if len(self.G2O) == 0 or len(self.O2G) == 0:
            raise IOError("WARNING: SNP pair file is EMPTY")

        return self.G2O, self.O2G, self.case_GWAS_set, self.case_outside_set, self.control_GWAS_set, self.control_outside_set

    def get_all_filenames(self):
        """
        read through filename file and create lists of case and control filenames
        """
        logging.debug(os.getcwd())
        f = open(self.case_file_file, "r")

        self.case_filenames = [self.filepath + line.strip() for line in f]
        f.close()

        f = open(self.control_file_file, "r")
        self.control_filenames = [self.filepath + line.strip() for line in f]
        f.close()

    @timeit
    def G_O_dict_maker(self):
        """
        creates case and control GWAS and outside genetic info data structures
        reads through all case and control files, updating data structures after every pass
        """
        self.case_GWAS = {}
        self.case_outside = {}
        self.control_GWAS = {}
        self.control_outside = {}

        self.get_all_filenames()

        skipped_files = []

        case_found_GWAS = set()
        case_found_outside = set()
        for case_file in self.case_filenames:
            if case_file in self.file_set:
                new_case_GWAS, new_case_outside, new_found_GWAS, new_found_outside = self.G_O_dict_maker_single(
                    case_file, self.case_GWAS_set, self.case_outside_set)
                self.case_GWAS.update(new_case_GWAS)
                self.case_outside.update(new_case_outside)
                case_found_GWAS.update(new_found_GWAS)
                case_found_outside.update(new_found_outside)
            else:
                skipped_files.append(case_file)

        control_found_GWAS = set()
        control_found_outside = set()
        for control_file in self.control_filenames:
            if control_file in self.file_set:
                new_control_GWAS, new_control_outside, new_found_GWAS, new_found_outside = self.G_O_dict_maker_single(
                    control_file, self.control_GWAS_set, self.control_outside_set)
                self.control_GWAS.update(new_control_GWAS)
                self.control_outside.update(new_control_outside)
                control_found_GWAS.update(new_found_GWAS)
                control_found_outside.update(new_found_outside)
            else:
                skipped_files.append(control_file)

        if len(skipped_files) > 0:
            logging.debug(
                "WARNING: The following files were skipped because they do not appear in the pairing file.\n" + str(skipped_files))

        error = False
        logging.debug(str(len(case_found_GWAS)) + " out of " +
                      str(len(self.case_GWAS_set)) + " GWAS rsIDs found in the case files.")
        not_found_GWAS = self.case_GWAS_set.difference(case_found_GWAS)
        if len(not_found_GWAS) > 0:
            logging.debug(
                "WARNING: The following GWAS rsIDs are found in the pairing file, but not in their case genotype file.\n" + str(not_found_GWAS))
            if not_found_GWAS == self.case_GWAS_set:
                error = True

        logging.debug(str(len(control_found_GWAS)) + " out of " +
                      str(len(self.control_GWAS_set)) + " GWAS rsIDs found in the control files.")
        not_found_GWAS = self.control_GWAS_set.difference(control_found_GWAS)
        if len(not_found_GWAS) > 0:
            logging.debug(
                "WARNING: The following GWAS rsIDs are found in the pairing file, but not in their control genotype file.\n" + str(not_found_GWAS))
            if not_found_GWAS == self.control_GWAS_set:
                error = True

        logging.debug(str(len(case_found_outside)) + " out of " +
                      str(len(self.case_outside_set)) + " outside rsIDs found in the case files.")
        not_found_outside = self.case_outside_set.difference(
            case_found_outside)
        if len(not_found_outside) > 0:
            logging.debug(
                "WARNING: The following outside rsIDs are found in the pairing file, but in their case genotype file.\n" + str(not_found_outside))
            if not_found_outside == self.case_outside_set:
                error = True

        logging.debug(str(len(control_found_outside)) + " out of " +
                      str(len(self.control_outside_set)) + " outside rsIDs found in the contol files.")
        not_found_outside = self.control_outside_set.difference(
            control_found_outside)
        if len(not_found_outside) > 0:
            logging.debug(
                "WARNING: The following outside rsIDs are found in the pairing file, but in their control genotype file.\n" + str(not_found_outside))
            if not_found_outside == self.control_outside_set:
                error = True

        if error:
            raise LookupError(
                "Check messages above. None of the GWAS/outside rsIDs were found in the .gen files.")

        return self.case_GWAS, self.case_outside, self.control_GWAS, self.control_outside

    def G_O_dict_maker_single(self, filename, GWAS_set, outside_set):
        """
        creates intermediate data structures for every run through a single file
        """
        file = open(filename, "r")

        found_GWAS = set()
        found_outside = set()
        GWAS_map = {}
        outside_map = {}

        line_num = 1
        for line in file:
            rsid = self.column_getter(line, self.rsid_col)
            check_key = (filename, rsid)
            if check_key in GWAS_set or check_key in outside_set:
                pos_1 = self.column_getter(line, self.snp_col)
                if "/" in pos_1:
                    pos_2 = pos_1[2]
                    pos_1 = pos_1[0]
                else:
                    pos_2 = self.column_getter(line, self.snp_col + 1)
                if len(pos_1) > 1 or len(pos_2) > 1:
                    continue

                poss_geno = possible_genotypes(pos_1, pos_2)
                head, alleles = self.column_cutter(line, self.col_cut)
                allele_list = alleles.split(self.delim)
                if self.is_triples:
                    if len(allele_list) % 3 != 0:
                        raise ValueError(
                            "Missing data. Incomplete genotype triple found in line " + str(line_num) + ".")
                    else:
                        allele_list = self.triplicate_converter(
                            allele_list, pos_1, pos_2)

                allele_list = list(map(self.clean_item_cached, allele_list))
                total = len(allele_list)
                if check_key in GWAS_set:
                    found_GWAS.add(check_key)
                    new_map = defaultdict(list)
                    new_map["alleles"] = poss_geno
                    new_map["total"] = total
                    for i in range(len(allele_list)):
                        curr_allele = allele_list[i]
                        if curr_allele != self.NA:
                            new_map[curr_allele].append(i)
                        else:
                            new_map["NA_indices"].append(i)
                    GWAS_map[self.swap_key_dict[check_key]] = new_map
                else:
                    found_outside.add(check_key)
                    outside_map[self.swap_key_dict[
                        check_key]] = allele_list, poss_geno
            line_num += 1
        file.close()

        return GWAS_map, outside_map, found_GWAS, found_outside

    def random_file_maker(self):
        Trans_Pipeline.random_file_maker(self)

    def random_file_maker_partial(self, outside_loc, GWAS_loc, genotypes):
        Trans_Pipeline.random_file_maker_partial(
            self, outside_loc, GWAS_loc, genotypes)

    def p_file_maker(self):
        Trans_Pipeline.p_file_maker(self)

    def p_file_maker_partial(self, outside_loc, GWAS_loc, genotypes):
        Trans_Pipeline.p_file_maker_partial(
            self, outside_loc, GWAS_loc, genotypes)