import pandas as pd
import ast
import os
import numpy as np
from scipy.spatial import distance
import sigfig



def get_alpaca_format_tree(simulations_path, tumour_id, normal_clone):
    tree_sim_format = pd.read_table(f"{simulations_path}/{tumour_id}_treeedges.tsv")[
        ["Parent", "Child"]
    ]
    tree_sim_format["Parent"] = "clone" + tree_sim_format["Parent"].astype(str)
    tree_sim_format["Child"] = "clone" + tree_sim_format["Child"].astype(str)
    # remove normal:
    tree_sim_format = tree_sim_format[tree_sim_format["Parent"] != normal_clone]
    MRCA = [
        x
        for x in tree_sim_format["Parent"].unique()
        if x not in tree_sim_format["Child"].unique()
    ]
    assert len(MRCA) == 1, "More than one MRCA in the tree"
    MRCA = MRCA[0]

    leaves = [
        x
        for x in tree_sim_format["Child"].unique()
        if x not in tree_sim_format["Parent"].unique()
    ]
    tree = []
    for branch in leaves:
        path = []
        current_clone = branch
        path.append(current_clone)
        while path[-1] != MRCA:
            current_parent = tree_sim_format[tree_sim_format["Child"] == current_clone][
                "Parent"
            ].values[0]
            path.append(current_parent)
            current_clone = current_parent
        path.reverse()
        tree.append(path)
    return tree


def compare_copynumber_profiles(TRUE_copynumbers,model_df,metric='hamming',copy_number_columns=['pred_CN_A','pred_CN_B']):
    model_df.sort_values(by='segment', inplace=True)
    TRUE_copynumbers.sort_values(by='segment', inplace=True)
    distance_func = getattr(distance, metric)
    results_dict = {}
    for true_clone in TRUE_copynumbers.clone.unique():
        true_clone_df = TRUE_copynumbers[TRUE_copynumbers.clone==true_clone]
        true_clone_unstacked = true_clone_df[copy_number_columns].unstack().values
        distances = {}
        # find most similar model clone:
        for model_clone in model_df.clone.unique():
            model_clone_df = model_df[model_df.clone==model_clone]
            assert sum((model_clone_df.segment.values) != (true_clone_df.segment.values)) == 0, 'Segments do not match'
            model_clone_unstacked = model_clone_df[copy_number_columns].unstack().values
            diff = distance_func(true_clone_unstacked, model_clone_unstacked)
            distances[model_clone] = diff
        closest_clone = min(distances, key=distances.get)
        closest_clone_distance = distances[closest_clone]
        results_dict[true_clone] = {
            "closest_clone": closest_clone,
            "closest_clone_distance": closest_clone_distance,
        }
    results_df = pd.DataFrame(results_dict).T.reset_index()
    return results_df


def set_mount_point():
    if os.getcwd().startswith("/Users"):
        mount_point_working = "/Volumes/lab-swantonc"
        mount_point_nemo = "/Volumes/proj-tracerx-lung"
    else:
        mount_point_working = "/nemo/lab/swantonc"
        mount_point_nemo = "/nemo/project/proj-tracerx-lung"
    return mount_point_working, mount_point_nemo


def cast_MEDICC_to_ALPACA(df):
    df["chrom"] = df["chrom"].str.replace("chr", "")
    df["segment"] = (
        df["chrom"] + "_" + df["start"].astype(str) + "_" + df["end"].astype(str)
    )
    df.rename(
        columns={
            "sample_id": "clone",
            "cn_a": "pred_CN_A",
            "cn_b": "pred_CN_B",
            "chrom": "chr",
            "start": "bin_start",
            "end": "bin_end",
        },
        inplace=True,
    )
    return df[
        ["clone", "pred_CN_A", "pred_CN_B", "segment", "chr", "bin_start", "bin_end"]
    ]


def get_HATCHET_proportions(HATCHET_output):
    HATCHET_proportions = HATCHET_output[
        [x for x in HATCHET_output.columns if "u_" in x]
    ].drop_duplicates()
    HATCHET_proportions = HATCHET_proportions.drop(columns="u_normal")
    HATCHET_proportions.columns = [x.split("_")[1] for x in HATCHET_proportions.columns]
    HATCHET_proportions = HATCHET_proportions.transpose()
    HATCHET_proportions = HATCHET_proportions.mean(axis=1)
    HATCHET_proportions_sum = HATCHET_proportions.sum()
    HATCHET_proportions_normalized = HATCHET_proportions / HATCHET_proportions_sum
    HATCHET_proportions_normalized = (
        HATCHET_proportions_normalized.reset_index().rename(
            columns={"index": "clone", 0: "proportion"}
        )
    )
    return HATCHET_proportions_normalized


def cast_HATCHET_to_ALPACA(HATCHET_output):
    HATCHET_output["segment"] = (
        HATCHET_output["#CHR"].astype(str)
        + "_"
        + HATCHET_output["START"].astype(str)
        + "_"
        + HATCHET_output["END"].astype(str)
    )
    HATCHET_output.drop(columns=["cn_normal"], inplace=True)
    HATCHET_output = HATCHET_output[
        [x for x in HATCHET_output.columns if "cn_" in x] + ["segment"]
    ].drop_duplicates()
    HATCHET_output.set_index("segment", inplace=True)
    HATCHET_alpaca_format = []
    for clone in HATCHET_output.columns:
        H_clone = HATCHET_output[[clone]].copy()
        H_clone["pred_CN_A"] = H_clone[clone].apply(lambda x: int(x.split("|")[0]))
        H_clone["pred_CN_B"] = H_clone[clone].apply(lambda x: int(x.split("|")[1]))
        H_clone.drop(columns=clone, inplace=True)
        H_clone["clone"] = clone.replace("cn_", "")
        HATCHET_alpaca_format.append(H_clone)
    HATCHET_alpaca_format = pd.concat(HATCHET_alpaca_format)
    HATCHET_alpaca_format.reset_index(inplace=True)
    return HATCHET_alpaca_format


def cast_TRUE_copynumbers_to_ALPACA(TRUE_copynumbers):
    TRUE_copynumbers = TRUE_copynumbers.drop(columns=["chr", "start", "end"])
    A_cols = [x for x in TRUE_copynumbers.columns if "A" in x]
    B_cols = [x for x in TRUE_copynumbers.columns if "B" in x]
    A_df_wide = TRUE_copynumbers[["segment"] + A_cols].copy()
    B_df_wide = TRUE_copynumbers[["segment"] + B_cols].copy()
    A_df_long = A_df_wide.melt(
        id_vars="segment", var_name="clone", value_name="pred_CN_A"
    )
    B_df_long = B_df_wide.melt(
        id_vars="segment", var_name="clone", value_name="pred_CN_B"
    )
    A_df_long["clone"] = A_df_long["clone"].str.replace("_A_true", "")
    B_df_long["clone"] = B_df_long["clone"].str.replace("_B_true", "")
    TRUE_copynumbers = pd.merge(A_df_long, B_df_long, on=["segment", "clone"])
    return TRUE_copynumbers


def get_proportions_whole_tumour(cp_table):
    cp_table.set_index("clone", inplace=True)
    clone_proportions = (
        cp_table.mean(axis=1).reset_index().rename(columns={0: "proportion"})
    )
    return clone_proportions


def get_ALPACA_cn_states(segment, ALPACA_output):
    ALPACA_CN_A = ALPACA_output[ALPACA_output["segment"] == segment]["pred_CN_A"].values
    ALPACA_CN_B = ALPACA_output[ALPACA_output["segment"] == segment]["pred_CN_B"].values
    return ALPACA_CN_A, ALPACA_CN_B


def get_HATCHET_cn_states(segment, HATCHET_output):
    H_segments = HATCHET_output[HATCHET_output["segment"] == segment]
    cols_to_drop = [x for x in H_segments.columns if "u_" in x] + ["SAMPLE"]
    H_segments = H_segments.drop(columns=cols_to_drop)
    assert (
        H_segments.duplicated().any()
    ), "Output has different clone copy numbers for the same segment"
    H_segments = H_segments.drop_duplicates()
    HATCHET_CN_cols = H_segments[[x for x in H_segments.columns if "cn_clone" in x]]
    HATCHET_CN_A = [int(x.split("|")[0]) for x in HATCHET_CN_cols.values[0]]
    HATCHET_CN_B = [int(x.split("|")[0]) for x in HATCHET_CN_cols.values[0]]
    return HATCHET_CN_A, HATCHET_CN_B


class CloneNameMismatchError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


def get_true_copy_numbers_per_POS(events_input_path, tree):
    unique_clones = set(sum(tree, []))
    events_input = pd.read_table(events_input_path).drop(columns="VAR_POS")
    clone_columns = [x for x in events_input.columns if x.startswith("(")]
    clone_columns_new_names = [f"clone{ast.literal_eval(x)[1]}" for x in clone_columns]
    events_input = events_input.rename(
        columns=dict(zip(clone_columns, clone_columns_new_names))
    )
    if set(unique_clones) != set(clone_columns_new_names):
        raise CloneNameMismatchError("Clone names in events table and tree don't match")
    for uc in unique_clones:
        if uc not in events_input.columns:
            events_input[uc] = np.nan
    clone_columns = [x for x in events_input.columns if "clone" in x]
    events_input_A = events_input.copy()
    events_input_B = events_input.copy()
    events_input_A[clone_columns] = events_input_A[clone_columns].applymap(
        lambda x: translateDescriptionToNumbers(x, "A")
    )
    events_input_B[clone_columns] = events_input_B[clone_columns].applymap(
        lambda x: translateDescriptionToNumbers(x, "B")
    )
    events_output_A = events_input_A.copy()
    events_output_B = events_input_B.copy()
    for clone in unique_clones:
        events_output_A.loc[:, clone] = getCloneCns(
            events_input_A, clone, tree, events_input_B
        )
        events_output_B.loc[:, clone] = getCloneCns(
            events_input_B, clone, tree, events_input_A
        )
    events_output_A = events_output_A.rename(
        columns={
            name_old: f"{name_old}_A_true"
            for name_old in events_output_A.columns
            if "clone" in name_old
        }
    )  # .drop(columns='diploid')
    events_output_B = events_output_B.rename(
        columns={
            name_old: f"{name_old}_B_true"
            for name_old in events_output_B.columns
            if "clone" in name_old
        }
    )  # .drop(columns='diploid')
    events_output = pd.merge(events_output_A, events_output_B)
    for clone_col in [c for c in events_output.columns if "clone" in c]:
        events_output.loc[(events_output[clone_col] < 0), clone_col] = 0
    return events_output


def translateDescriptionToNumbers(listString, allele):
    # other_allele = 'A' if allele=='B' else 'B'
    if listString != listString:
        return 0
    else:
        components = ast.literal_eval(listString)
        description = {
            "A": {
                "cnloh_A": "*0",
                "cnloh_B": "+x",
                "snv": "+0",
                "wgd": "*2",
                "loss_A": "-1",
                "gain_both": "+1",
                "gain_A": "+1",
            },
            "B": {
                "cnloh_A": "+x",
                "cnloh_B": "*0",
                "snv": "+0",
                "wgd": "*2",
                "loss_B": "-1",
                "gain_both": "+1",
                "gain_B": "+1",
            },
        }
        description_alleles = description[allele]
        components = [
            description_alleles[v]
            for v in components
            if v in list(description_alleles.keys())
        ]
        """if '+x' in components_this_allele:
            description_alleles_other = description[other_allele]
            components_other_allele = [description_alleles_other[v] for v in components if v in list(description_alleles_this.keys())]
        """
        return str(components)


def getCloneCns(allele_df, clone, tree, other_allele_df):
    ancestors = getAncestors(clone, tree)
    clone_sub_df = allele_df[[x for x in ancestors if x in allele_df.columns]].copy()
    if "+x" in str(clone_sub_df.unstack().to_list()):
        clone_sub_df_with_x = clone_sub_df[
            clone_sub_df.apply(lambda x: x.str.contains("\+x", regex=True))
            .fillna(False)
            .sum(axis=1)
            .astype(bool)
        ]
        clone_sub_df_other = other_allele_df[
            [x for x in ancestors if x in other_allele_df.columns]
        ]
        clone_sub_df_other_where_x = clone_sub_df_other[
            clone_sub_df.apply(lambda x: x.str.contains("\+x", regex=True))
            .fillna(False)
            .sum(axis=1)
            .astype(bool)
        ]
        clone_cns_other = clone_sub_df_other_where_x.apply(
            lambda x: evaluateStringExpressionsForX(x), axis=1
        )
        for i, (ind, value) in enumerate(clone_cns_other.iteritems()):
            clone_sub_df.iloc[ind] = clone_sub_df.iloc[ind].apply(
                lambda x: replaceX(x, value)
            )
    clone_cns = clone_sub_df.apply(
        lambda x: evaluateStringExpressions(x), axis=1
    )  # clone_sub_df.sum(axis=1)
    return clone_cns


def evaluateStringExpressions(row):
    vals = flat_list([ast.literal_eval(x) for x in row if x != 0])
    final_value = "1"
    for v in vals:
        expression = final_value + v
        final_value = str(eval(expression))
    return int(final_value)


def getAncestors(clone, tree):
    branch = [b for b in tree if clone in b][0]
    ancestors = branch[: branch.index(clone) + 1]
    return ancestors


def evaluateStringExpressionsForX(row):
    vals = flat_list([ast.literal_eval(x) for x in row if x != 0])
    vals = vals[: vals.index("*0")]
    final_value = "1"
    for v in vals:
        expression = final_value + v
        final_value = str(eval(expression))
    return int(final_value)


def flat_list(target_list):
    if (isinstance(target_list, list)) & (len(target_list) > 0):
        if isinstance(target_list[0], list):
            return [item for sublist in target_list for item in sublist]
    else:
        return target_list


def replaceX(cell, v):
    if isinstance(cell, str):
        if "+x" in cell:
            if "+x" in cell:
                cell = cell.replace("x", str(v))
    return cell


# make true copy numbers df in ALPACA format
def make_true_cn_alpaca_format(
    HATCHET_output, TRUE_copynumbers, events_output, tumour_id
):
    seg_res = []
    diploid_segments = []
    tetraploid_segments = []
    for s in HATCHET_output["segment"].unique():
        chromosome = int(s.split("_")[0])
        start = int(s.split("_")[1])
        end = int(s.split("_")[2])
        TRUE_segment_bins = TRUE_copynumbers[
            (TRUE_copynumbers["CHR"] == chromosome)
            & (TRUE_copynumbers["START"] >= start)
            & (TRUE_copynumbers["END"] <= end)
        ]
        TRUE_segment_bins_nonans = TRUE_segment_bins.dropna(subset=["POS"])
        if TRUE_segment_bins_nonans.shape[0] == 0:
            if (TRUE_segment_bins.TRUE_TOTALCN_TUMOUR == 2).all():
                diploid_segments.append(s)
            elif (TRUE_segment_bins.TRUE_TOTALCN_TUMOUR == 4).all():
                tetraploid_segments.append(s)
            continue
        POS = TRUE_segment_bins_nonans.POS.unique()[0]
        segment_true = events_output[events_output.POS.isin([POS])].iloc[0]
        clone_names = [col for col in segment_true.index if "clone" in col]
        clone_values_A = [segment_true[col] for col in clone_names if "_A_true" in col]
        clone_values_B = [segment_true[col] for col in clone_names if "_B_true" in col]
        clone_names = [col.split("_")[0] for col in clone_names if "_A_true" in col]
        segment_true_ALPACA = pd.DataFrame(
            {
                "clone": clone_names,
                "pred_CN_A": clone_values_A,
                "pred_CN_B": clone_values_B,
            }
        )
        segment_true_ALPACA["segment"] = s
        segment_true_ALPACA["tumour_id"] = tumour_id
        seg_res.append(segment_true_ALPACA)
        normal_template = segment_true_ALPACA.copy()
        for diploid_segment in diploid_segments:
            this_segment = normal_template.copy()
            this_segment["segment"] = diploid_segment
            this_segment["pred_CN_A"] = 1
            this_segment["pred_CN_B"] = 1
            seg_res.append(this_segment)
        for tetraploid_segment in tetraploid_segments:
            this_segment = normal_template.copy()
            this_segment["segment"] = tetraploid_segment
            this_segment["pred_CN_A"] = 2
            this_segment["pred_CN_B"] = 2
            seg_res.append(this_segment)
    true_cn_ALPACA = pd.concat(seg_res)
    true_cn_ALPACA = true_cn_ALPACA.drop_duplicates().reset_index(drop=True)
    # true_cn_ALPACA.to_csv(f'/camp/home/pawlikp/CN-CCF/publication/output/H2/default/patient_outputs/{tumour_id}/true_copynumbers_true_clones.csv', index=False)
    return true_cn_ALPACA


# old version, proportions only:
def score_segment_with_proportions_OLD(
    segment, TRUE_copynumbers, true_proportions, model_output, model_proportions
):
    # for each segment and each true copy number state, compare the clone proportions of clones with such state in true
    # and in model output. If they are exactly the same, the error is 0.
    seg_TRUE_copynumbers = TRUE_copynumbers[TRUE_copynumbers.segment == segment][
        ["clone", "pred_CN_A", "pred_CN_B"]
    ]
    seg_model_copynumbers = model_output[model_output.segment == segment][
        ["clone", "pred_CN_A", "pred_CN_B"]
    ]
    seg_TRUE_copynumbers = seg_TRUE_copynumbers.merge(true_proportions, on="clone")
    seg_model_copynumbers = seg_model_copynumbers.merge(model_proportions, on="clone")
    proportions_by_state_A_true = (
        seg_TRUE_copynumbers[["proportion", "pred_CN_A"]].groupby("pred_CN_A").sum()
    )
    proportions_by_state_B_true = (
        seg_TRUE_copynumbers[["proportion", "pred_CN_B"]].groupby("pred_CN_B").sum()
    )
    proportions_by_state_A_model = (
        seg_model_copynumbers[["proportion", "pred_CN_A"]].groupby("pred_CN_A").sum()
    )
    proportions_by_state_B_model = (
        seg_model_copynumbers[["proportion", "pred_CN_B"]].groupby("pred_CN_B").sum()
    )
    proportions_by_state_A = proportions_by_state_A_true.merge(
        proportions_by_state_A_model
    )
    proportions_by_state_A = proportions_by_state_A_true.merge(
        proportions_by_state_A_model,
        left_index=True,
        right_index=True,
        how="left",
        suffixes=("_true", "_model"),
    ).fillna(0)
    proportions_by_state_A["error"] = abs(
        proportions_by_state_A["proportion_true"]
        - proportions_by_state_A["proportion_model"]
    )
    proportions_by_state_B = proportions_by_state_B_true.merge(
        proportions_by_state_B_model,
        left_index=True,
        right_index=True,
        how="left",
        suffixes=("_true", "_model"),
    ).fillna(0)
    proportions_by_state_B["error"] = abs(
        proportions_by_state_B["proportion_true"]
        - proportions_by_state_B["proportion_model"]
    )
    mean_error_A = proportions_by_state_A["error"].mean()
    mean_error_B = proportions_by_state_B["error"].mean()
    mean_error = (mean_error_A + mean_error_B) / 2
    return mean_error


# new version, with wasserstein distance:
def score_segment_with_proportions(
    TRUE_seg,
    true_proportions,
    ALPACA_seg,
    conipher_proportions,
    HATCHET_seg,
    HATCHET_proportions,
):
    allele_res = []
    for allele in ["A", "B"]:
        proportions_by_state_ALPACA = (
            ALPACA_seg[["clone", f"pred_CN_{allele}", "segment"]]
            .merge(conipher_proportions)[[f"pred_CN_{allele}", "proportion"]]
            .groupby(f"pred_CN_{allele}")
            .sum()
        )
        proportions_by_state_ALPACA.reset_index(inplace=True)
        proportions_by_state_ALPACA.rename(
            columns={f"pred_CN_{allele}": "cn"}, inplace=True
        )
        proportions_by_state_TRUE = (
            TRUE_seg[["clone", f"pred_CN_{allele}", "segment"]]
            .merge(true_proportions)[[f"pred_CN_{allele}", "proportion"]]
            .groupby(f"pred_CN_{allele}")
            .sum()
        )
        proportions_by_state_TRUE.reset_index(inplace=True)
        proportions_by_state_TRUE.rename(
            columns={f"pred_CN_{allele}": "cn"}, inplace=True
        )
        proportions_by_state_HATCHET = (
            HATCHET_seg[["clone", f"pred_CN_{allele}", "segment"]]
            .merge(HATCHET_proportions)[[f"pred_CN_{allele}", "proportion"]]
            .groupby(f"pred_CN_{allele}")
            .sum()
        )
        proportions_by_state_HATCHET.reset_index(inplace=True)
        proportions_by_state_HATCHET.rename(
            columns={f"pred_CN_{allele}": "cn", "proportion": "proportion_hatchet"},
            inplace=True,
        )
        proportions_by_state = proportions_by_state_TRUE.merge(
            proportions_by_state_ALPACA,
            on="cn",
            how="left",
            suffixes=("_true", "_alpaca"),
        ).fillna(0)
        proportions_by_state = proportions_by_state.merge(
            proportions_by_state_HATCHET, on="cn", how="left"
        ).fillna(0)
        proportions_by_state["error_ALPACA"] = abs(
            proportions_by_state["proportion_true"]
            - proportions_by_state["proportion_alpaca"]
        )
        proportions_by_state["error_HATCHET"] = abs(
            proportions_by_state["proportion_true"]
            - proportions_by_state["proportion_hatchet"]
        )
        proportions_by_state["allele"] = allele
        allele_res.append(proportions_by_state)
    return pd.concat(allele_res)


def is_multiple_of(list1, list2):
    if len(list1) > len(list2):
        list1, list2 = list2, list1
    if len(list2) % len(list1) != 0:
        return False
    repeat_count = len(list2) // len(list1)
    return list1 * repeat_count == list2


def filter_failed_solutions(df, expected_clones):
    filtered_dfs = []
    for c in df.complexity.unique():
        df_c = df[df.complexity == c]
        df = df[df.complexity != c]
        if list(df_c.clone) == expected_clones:
            filtered_dfs.append(df_c)
        # in some cases, multiple solutions have the same complexity:
        elif is_multiple_of(expected_clones, list(df_c.clone)):
            solutions_count_in_df_c = len(df_c.clone) // len(expected_clones)
            for i in range(solutions_count_in_df_c):
                filtered_dfs.append(
                    df_c.iloc[i * len(expected_clones) : (i + 1) * len(expected_clones)]
                )
    return pd.concat(filtered_dfs)


def add_solution_identifier(combined, all_solutions_seg, tumour_id, segment):
    # solutions lack a clear identifier, so split the data frame into sub dataframes using the shape of the optimal solution:
    alpaca_solution_clones = combined.query(
        "tumour_id == @tumour_id & segment == @segment"
    )
    solution_count = all_solutions_seg.shape[0] / alpaca_solution_clones.shape[0]
    # check if solution count is an integer:
    assert solution_count.is_integer(), "Solution count is not an integer"
    solution_count = int(solution_count)
    solution_vector = (
        list(range(1, solution_count + 1)) * alpaca_solution_clones.shape[0]
    )
    solution_vector.sort()
    all_solutions_seg.loc[:, "solution"] = solution_vector
    return all_solutions_seg


def get_hd_hatchet_alpaca_colours(palette):
    import json

    with open(palette, "r") as f:
        palette = json.load(f)
    catgorical_palette = palette["categorical"]
    colours = {}
    colours["hd"] = catgorical_palette["c12"]
    colours["hatchet"] = catgorical_palette["c8"]
    colours["hatchet2"] = catgorical_palette["c5"]
    colours["alpaca"] = catgorical_palette["c2"]
    return colours


def get_true_wgd(tumour_id, simulations_input_path):
    tree_edges = pd.read_csv(f'{simulations_input_path}/{tumour_id}/sim/{tumour_id}_treeedges.tsv', sep='\t')
    true_wgd = {'clonal':[f'clone{x}' for x in tree_edges[(tree_edges['WGD']>0)&(tree_edges['Truncal']=='Yes')]['Child'].values], 'subclonal':[f'clone{x}' for x in tree_edges[(tree_edges['WGD']>0)&(tree_edges['Truncal']=='No')]['Child'].values]}
    return true_wgd


def get_M2_wgd(tumour_id, hatchet_output):
    tumour_dir = f'{hatchet_output}/{tumour_id}'
    medicc_dir = f'{tumour_dir}/medicc'
    cn_profiles = pd.read_csv(f'{medicc_dir}/{tumour_id}_final_cn_profiles.tsv', sep='\t')
    wgd_status_df = cn_profiles[['sample_id','is_wgd']].drop_duplicates()
    wgd_status_df = wgd_status_df[wgd_status_df['is_wgd']]
    m2_wgd = {'clonal':[], 'subclonal':[]}
    if wgd_status_df.shape[0] == 0:
        # lack of wgd encoded as empty list
        return m2_wgd
    else:
        wgd_clones = list(wgd_status_df['sample_id'].values)
        # sometimes M2 encodes the internal clone 1 as 'internal_1' and sometimes as 'internal1'
        if 'internal_1' in wgd_clones:
            wgd_clones = [x.replace('_','') for x in wgd_clones]
        if ('internal1' in wgd_clones):
            m2_wgd['clonal'].append('internal1')
        wgd_subclones = [x for x in wgd_clones if x not in m2_wgd['clonal']]
        m2_wgd['subclonal'] = wgd_subclones
    return m2_wgd


def get_ALPACA_wgd(tumour_id, alpaca_wgd_calls_df):
    tumour_wgd = alpaca_wgd_calls_df[alpaca_wgd_calls_df['tumour_id']==tumour_id]
    tumour_wgd_clones = tumour_wgd[tumour_wgd['GD']>0]
    alpaca_wgd = {'clonal': list(tumour_wgd_clones[tumour_wgd_clones.MRCA].clones.values), 'subclonal':list(tumour_wgd_clones[~tumour_wgd_clones.MRCA].clones.values)}
    return alpaca_wgd


def pval_to_annotations(pvalue):
    if pvalue < 0.001:
        formatted = sigfig.round(pvalue, sigfigs=3, notation='scientific')
        base, exponent = formatted.split('E')
        symbol = f"P = {base}×10<sup>{int(exponent)}</sup>"
    else:
        symbol = f"P = {sigfig.round(pvalue, sigfigs=3)}"
    return symbol

