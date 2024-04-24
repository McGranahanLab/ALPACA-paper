#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
from python_functions import get_ancestors, get_descendants, find_parent, read_tree_json, get_edge


def is_child_of_gd(tree, r):
    c = []
    gd_clones = r.clones[r.GD > 0]
    for clone in r['clones']:
        ancestors = [x for x in get_ancestors(clone, tree) if x != clone]
        if bool(sum([x in ancestors for x in gd_clones])):
            c.append('child')
        else:
            c.append('')
    r['isChildOfGD'] = c
    return r


def has_gd_child(tree, r):
    c = []
    gd_clones = r.clones[r.GD > 0]
    for clone in r['clones']:
        descendants = [x for x in get_descendants(clone, tree) if x != clone]
        if bool(sum([x in descendants for x in gd_clones])):
            c.append('parent')
        else:
            c.append('')
    r['isParentOfGD'] = c
    return r


def call_gd_on_edge(edge, grandparent, alpaca_output, MRCA, events_per_clone):
    """
    for each clone (child) calculate the rations between copy number states of:
    child and parent
    child and grandparent
    if this ratio is > 1.5, the child is classified as GD
    Additionally, if the ratio between the child and parent is below 1.5 but the ratio between the child and grandparent is above 1.5,
    the child is classified as GD unless the parent is already classified as GD
    """
    parent, child = edge
    alpaca_output['clone'] = alpaca_output['clone'].astype(str)
    alpaca_output['seg_len'] = alpaca_output.segment.str.split('_', expand=True)[2].astype(int) - alpaca_output.segment.str.split('_', expand=True)[1].astype(int)
    parent_frame = alpaca_output[alpaca_output.clone == parent]
    child_frame = alpaca_output[alpaca_output.clone == child]
    grandparent_frame = alpaca_output[alpaca_output.clone == grandparent]
    e = pd.DataFrame()
    C = child_frame.pred_CN_A
    P = parent_frame.pred_CN_A
    G = grandparent_frame.pred_CN_A
    ratio_CP = C.values / P.values
    ratio_CG = C.values / G.values
    ratio_CP[np.isnan(ratio_CP)] = 1
    ratio_CG[np.isnan(ratio_CG)] = 1
    e['parent_child_ratio'] = [np.average(a=ratio_CP, weights=parent_frame.seg_len)]
    if parent != 'clone100':
        e['grandparent_child_ratio'] = [np.average(a=ratio_CG, weights=parent_frame.seg_len)]
    e['clones'] = [child]
    if child == MRCA:
        e['MRCA'] = True
    else:
        e['MRCA'] = False
    e['GD'] = False
    clone_is_gd = (e['parent_child_ratio'].round() > 1).values[0]
    e.loc[e['parent_child_ratio'].round() > 1, 'GD'] = True
    if (parent != 'clone100') and (len(events_per_clone) > 0) and (not clone_is_gd):
        parent_was_gd, grandparent_was_gd = False, False
        try:
            parent_was_gd = bool(events_per_clone.loc[events_per_clone.clones == parent, 'GD'].iloc[0])
        except IndexError:
            pass
        try:
            grandparent_was_gd = bool(events_per_clone.loc[events_per_clone.clones == grandparent, 'GD'].iloc[0])
        except IndexError:
            pass
        # this is to capture the cases with gradual increase, i.e. grandparent at 1, parent at 1.5, child at 2
        if (not parent_was_gd) & (not grandparent_was_gd):
            child_should_be_gd = ((e['parent_child_ratio'].round() < 2) & (e['grandparent_child_ratio'].round() > 1)).values[0]
            if child_should_be_gd:
                e['GD'] = True
                e['GD_move_to'] = parent
    
    e['ploidy_A'] = np.average(a=child_frame.pred_CN_A, weights=child_frame.seg_len)
    e['ploidy_B'] = np.average(a=child_frame.pred_CN_B, weights=child_frame.seg_len)
    e['ploidy'] = np.average(a=child_frame.pred_CN_A + child_frame.pred_CN_B, weights=child_frame.seg_len)
    
    e['LOH_fraction_A'] = sum((child_frame.pred_CN_A == 0) * child_frame.seg_len) / sum(child_frame.seg_len)
    e['LOH_fraction_B'] = sum((child_frame.pred_CN_B == 0) * child_frame.seg_len) / sum(child_frame.seg_len)
    e['LOH_fraction'] = e['LOH_fraction_A'] + e['LOH_fraction_B']
    return e

    
    
def call_wgd(alpaca_output:pd.DataFrame, tree:list):
    alpaca_output = alpaca_output[['tumour_id', 'clone', 'pred_CN_A', 'pred_CN_B', 'segment']].drop_duplicates()
    MRCA = tree[0][0]
    if MRCA != 'diploid':
        # add 'diploid' to the tree
        for i, b in enumerate(tree):
            tree[i] = ['diploid'] + b
    edges, _ = get_edge(tree)
    events_per_clone = pd.DataFrame()
    if not 'diploid' in alpaca_output.clone.unique():
        # create diploid dataframe
        diploid_frame = alpaca_output[alpaca_output.clone == MRCA].copy()
        diploid_frame['pred_CN_A'] = 1
        diploid_frame['pred_CN_B'] = 1
        diploid_frame['clone'] = 'diploid'
        # add diploid to the alpaca_output
        alpaca_output = pd.concat([alpaca_output, diploid_frame])
    alpaca_output['chr'] = alpaca_output.segment.str.split('_', expand=True)[0]
    alpaca_output['start'] = alpaca_output.segment.str.split('_', expand=True)[1].astype(int)
    for edge in edges:
        print(edge)
        if (edge[0] == MRCA) or (edge[0] == 'diploid'):
            grandparent = 'diploid'
        else:
            grandparent = find_parent(edge[0], tree)
        e = call_gd_on_edge(edge, grandparent, alpaca_output, MRCA, events_per_clone)
        events_per_clone = pd.concat([events_per_clone, e])
    r = events_per_clone
    r['GD'] = r['GD'].astype(int)
    if 'GD_move_to' in r.columns:
        clones_to_move_up_the_tree = r[~r['GD_move_to'].isna()].clones
        for clone in clones_to_move_up_the_tree:
            p = r.loc[r.clones == clone, 'GD_move_to'].values[0]
            r.loc[r.clones == clone, 'GD'] = 0
            r.loc[r.clones == p, 'GD'] = 1
    if r.GD.sum() > 1:
        second_gd_candidates = [x for x in r.clones[r.GD == 1] if len(set([y for y in get_ancestors(x, tree) if y != x]).intersection(set(list(r.clones[r.GD == 1])))) > 0]
        for c in second_gd_candidates:
            if r.loc[r.clones == c, 'parent_child_ratio'].iloc[0].round() > 1:
                if r.loc[r.clones == c, 'ploidy_A'].iloc[0] > 3.5:
                    r.loc[r.clones == c, 'GD'] = 2
                else:
                    r.loc[r.clones == c, 'GD'] = 0
    for branch in tree:
        # everything after second GD is classified as gd child:
        second_gd_clone = r.loc[r.GD == 2, 'clones'].values
        if len(second_gd_clone) > 0:
            second_gd_clone_descendants = get_descendants(r.loc[r.GD == 2, 'clones'].values[0], tree)
            branch_GD_clones = [x for x in r.loc[(r.clones.isin(branch)) & (r.GD > 0)].clones if x not in second_gd_clone_descendants + [r.loc[r.GD == 2, 'clones'].iloc[0]]]
            r.loc[(r.clones.isin(second_gd_clone_descendants) & (r.GD == 1)), 'GD'] = 0
            # if there are two sequential clones both classified as first GD, filter all but one
            if len(branch_GD_clones) > 1:
                max_ploidy = r.loc[r.clones.isin(branch) & r.GD == 1].ploidy_A.max()
                max_ploidy_clone = r.loc[(r.clones.isin(branch)) & (r.GD == 1) & (r.ploidy_A == max_ploidy)].clones.values[0]
                competing_clones = [x for x in branch if x != max_ploidy_clone]
                r.loc[(r.clones.isin(competing_clones)) & (r.GD == 1), 'GD'] = 0
    r['tumour_id'] = alpaca_output.tumour_id.unique()[0]
    r = is_child_of_gd(tree, r)
    r = has_gd_child(tree, r)
    return r


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_data_directory",type=str, required=True,help='Directory where input data is stored. Should contain subdirectories for each tumour')
    parser.add_argument('--cohort_results_file', type=str, required=True, help='Results for all tumours')
    args = parser.parse_args()

    input_data_directory = args.input_data_directory
    cohort_results_file = args.cohort_results_file
    cohort_results = pd.read_csv(cohort_results_file)
    wgd_per_tumour_outputs = []
    for tumour_id, df in cohort_results.groupby('tumour_id'):
        print(tumour_id)
        tumour_dir = f'{input_data_directory}/{tumour_id}'
        tree = read_tree_json(f'{tumour_dir}/tree_paths.json')
        wgd_per_tumour_outputs.append(call_wgd(df, tree))
    df_out = pd.concat(wgd_per_tumour_outputs)
    df_out.to_csv(f'wgd_calls.csv')


if __name__ == "__main__":
    main()