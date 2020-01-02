import pandas as pd
import numpy as np

PIPE_TYPE = None

class MtcTable(object):

    """Summary
    Return an Mtc_table
    Handles whether the table actually exist or not
    Attributes:
            mtc_table (TYPE): Description
            threshold_col (TYPE): Description
    """

    def __init__(self, mtc_file_path, threshold_col):
        self.mtc_table = self.__get_table__(mtc_file_path)
        self.threshold_col = threshold_col

    @staticmethod
    def __get_table__(path):
        mtc = pd.read_csv(path, index_col=0, sep = "\t")
        mtc = mtc.sort_index()
        return mtc

    def get_threshold(self, GWAS_rsid):
        return self.mtc_table.loc[GWAS_rsid, self.threshold_col]

    @staticmethod
    def process_df(dictionary, dict_name, common_keys):
        from pipelines import Trans_Case_Control_Pipeline
        if PIPE_TYPE == Trans_Case_Control_Pipeline:
            processed_dict = {
                list(zip(*k))[1]: v for k, v in dictionary.items()}
            df = pd.DataFrame.from_dict(processed_dict, orient="index")
        else:
            df = pd.DataFrame.from_dict(dictionary, orient="index")
            df = df.loc[list(common_keys), :]

        #df.to_csv("process_df_{}".format(dict_name), sep="\t")

        non_geno_cols = ['NA_indices', "total"]
        non_geno_df = df[non_geno_cols].reset_index()
        non_geno_df["NA_indices"] = non_geno_df["NA_indices"].apply(len)
        non_geno_df["{}_total_no_NA".format(dict_name)] = non_geno_df[
            "total"] - non_geno_df['NA_indices']
        non_geno_df.rename(columns=dict(zip(non_geno_cols, ["{}_{}".format(
            dict_name, col) for col in non_geno_cols])), inplace=True)

        df.drop(columns=non_geno_cols, inplace=True)
        df = df.stack().reset_index()
        df = pd.merge(df, non_geno_df, on=["level_0", "level_1"])

        array_col_name = "{}_array".format(dict_name)
        df.rename(columns=dict(zip(df.columns[0:4], [
                  "GWAS_rsid", "outside_rsid", "geno_combo", array_col_name])), inplace=True)

        df[["GWAS_geno", "outside_geno"]] = pd.DataFrame(
            df["geno_combo"].values.tolist(), index=df.index)
        df = df.sort_values(
            by=["GWAS_rsid", "outside_rsid", "GWAS_geno", "outside_geno"])
        df["{}_count".format(dict_name)] = df[array_col_name].apply(len)
        df.drop(columns="geno_combo", inplace=True)

        return df

    @staticmethod
    def create_mtc_table(df, total_case_unfiltered, total_control_unfiltered):

        found_pairs = df[["GWAS_rsid", "outside_rsid"]].drop_duplicates().shape[
            0]

        df = df.query("case_count!=0 and control_count!=0 ")
        df = df.query(
            "(case_count/case_total_no_NA) > 0.01 and (control_count/control_total_no_NA) >0.01")
        df = df.sort_values(
            by=["GWAS_rsid", "outside_rsid", "GWAS_geno", "outside_geno"])
        geno_counts = df.groupby(["GWAS_rsid", "outside_rsid"]).size(
        ).reset_index(name="non_zero_geno_combo_counts")

        filter_1 = pd.merge(df, geno_counts, how="inner").query("non_zero_geno_combo_counts > 4")

        # made these a parameter of [sample_out file length - 2] (2 lines of
        # header): total_case_unfiltered, total_control_unfiltered
        filtered_after_1 = geno_counts.query("not(non_zero_geno_combo_counts > 4)")

        filter_2 = filter_1.query("(case_total_no_NA/{}) > 0.1 & (control_total_no_NA/{}) > 0.1".format(
            total_case_unfiltered, total_control_unfiltered))

        filtered_after_2 = filter_1.query(f"not ((case_total_no_NA/{total_case_unfiltered}) > 0.1 & (control_total_no_NA/{total_control_unfiltered}) > 0.1)")

        # count all GWAS_Outside genotype combinations per GWAS_rsid (for each
        # GWAS rsid summing/marginalizing over all the geno combo of all
        # outside rsid)
        mtc_table_created = filter_2.groupby("GWAS_rsid").size().reset_index(
            name="num_geno_combo_per_GWAS_rsid_filtered_more_than_4")
        mtc_table_created["threshold"] = 0.05 / mtc_table_created.iloc[:, 1]
        return mtc_table_created, found_pairs, filter_1, filter_2, filtered_after_1, filtered_after_2

    @classmethod
    def make_mtc_table_from_dict(cls, case_combined_dict, control_combined_dict, total_case_unfiltered, total_control_unfiltered):
        common_keys = set(case_combined_dict.keys()) & set(
            control_combined_dict.keys())

        control_combined_df = cls.process_df(
            control_combined_dict, "control", common_keys)
        # make df with the following columns: ['GWAS_rsid', 'outside_rsid',
        # 'case_array', 'case_NA_indices', 'case_total', 'case_total_no_NA',
        # 'GWAS_geno', 'outside_geno', 'case_count']
        case_combined_df = cls.process_df(
            case_combined_dict, "case", common_keys)

        case_control_combined_df = case_combined_df.merge(
            control_combined_df, on=["GWAS_rsid", "outside_rsid", "GWAS_geno", "outside_geno"])

        case_control_combined_df.to_csv("case_control_combined_df", sep="\t")

        mtc_table, found_pairs, filter_1, filter_2, filtered_after_1, filtered_after_2 = cls.create_mtc_table(
            case_control_combined_df, total_case_unfiltered, total_control_unfiltered)

        return mtc_table, found_pairs, filter_1, filter_2, filtered_after_1, filtered_after_2


class StepwiseFilter(object):
    z_dict = {99: 2.576, 98: 2.326, 95: 1.96, 90: 1.645, 80: 1.28}

    def __init__(self, mtc_table, z_threshold):
        self.mtc_table = mtc_table
        self.z = self.z_dict[int(z_threshold)]

    def get_c_i(self, iterations, p_value):

        factor = self.z * np.sqrt((p_value * (1 - p_value)) / iterations)
        c_i_neg = p_value - factor
        c_i_pos = p_value + factor

        return c_i_neg, c_i_pos

    def get_compare_info_list(self, iterations, p_value, GWAS_rsid, outside_cause_higher_or_lower_risk):

        # in other pipelines, GWAS_id is a string, in TRANS_CC pipeline, it is
        # a tuple of the form:
        # ('ALL_MS_impute2_chr16.gen,ALL_controls_58C_NBS_WTC2_impute2_chr16.gen',
        # 'rs13333054')
        from pipelines import Trans_Case_Control_Pipeline
        if PIPE_TYPE == Trans_Case_Control_Pipeline:
            GWAS_rsid = GWAS_rsid[1]
        # logging.debug(PIPE_TYPE)
        #logging.debug("GWAS_rsid: {}".format(GWAS_rsid))

        # handles p values that are zeros
        p_value_no_zero = p_value if p_value != 0 else float(
            1) / float(iterations)
        c_i_neg, c_i_pos = self.get_c_i(iterations, p_value_no_zero)
        mtc_threshold = self.mtc_table.get_threshold(GWAS_rsid)
        need_more_perm = (mtc_threshold >= c_i_neg) & (
            mtc_threshold <= c_i_pos)

        if not need_more_perm:
            if mtc_threshold < c_i_neg:
                status = "non_sig"
            else:
                status = "sig_{}".format(outside_cause_higher_or_lower_risk)
        else:
            status = "more_perm"

        return p_value_no_zero, c_i_neg, c_i_pos, mtc_threshold, need_more_perm, status