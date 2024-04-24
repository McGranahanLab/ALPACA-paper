#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
from python_functions import get_chr_table, flat_list


def add_abs_coordinates(df,chr_table_path):
    chr_table = get_chr_table(chr_table_path)
    df = df.merge(chr_table, on='chr')
    df['abs_start'] = df['start'] + df['shift']
    df['abs_end'] = df['end'] + df['shift']
    df.drop(columns=['shift','len','cumsum','ticks'], inplace=True)
    return df


def find_consensus_segmentation(list_of_dfs):
    # input: list of dataFrames with segment_name, abs_start and abs_end fields
    # output: list of dataFrames with consensus_segment, consensus_abs_start, consensus_abs_end fields + segment to merge with original dataframe
    # test data:
    # df1 = pd.DataFrame({'segment': ['0_100'], 'abs_start': [0], 'abs_end': [100], 'meta': ['x']})
    # df2 = pd.DataFrame({'segment': ['0_50', '51_100'], 'abs_start': [0, 51], 'abs_end': [50, 100], 'meta': ['x','y']})
    # df3 = pd.DataFrame({'segment': ['0_30', '31_70', '71_100'], 'abs_start': [0, 31, 71], 'abs_end': [30, 70, 100],'meta': ['x','y','z']})
    # list_of_dfs = [df1, df2, df3]
    # check if required columns are in all input dfs:
    for c in ['abs_start', 'abs_end', 'segment']:
        for d in list_of_dfs:
            assert (c in d.columns)
    # out
    output_dfs = []
    # breakpoint is defined as 'ends', excluding the last end:
    new_breakpoints = list(set(flat_list([(list(x['abs_start']) + list(x['abs_end'])) for x in list_of_dfs])))
    new_breakpoints = sorted(new_breakpoints)
    for i, d in enumerate(list_of_dfs):
        d = d.sort_values('abs_start', ascending=True)
        new_starts = new_breakpoints[:-1]
        new_ends = new_breakpoints[1:]
        new_segments = pd.DataFrame({'consensus_abs_start': new_starts, 'consensus_abs_end': new_ends})
        new_segments['consensus_segment'] = new_segments['consensus_abs_start'].astype(str) + '_' + new_segments['consensus_abs_end'].astype(str)
        new_segments['len'] = new_segments['consensus_abs_end'] - new_segments['consensus_abs_start']
        new_segments = new_segments[new_segments['len'] > 5]
        new_segments['segment'] = new_segments.apply(lambda row: find_original_matching_segment(row, d), axis=1)
        new_segments = pd.merge(new_segments, d, on='segment')
        output_dfs.append(new_segments)
    list_of_all_segments = [set(x['consensus_segment']) for x in output_dfs]
    common_segments = list(list_of_all_segments.pop().intersection(*list_of_all_segments))
    output_dfs = [df[df['consensus_segment'].isin(common_segments)] for df in output_dfs]
    # all output dataframes should now have the same number of segments
    assert (len(set([len(x) for x in output_dfs])) == 1)
    return output_dfs


def find_original_matching_segment(row, df):
    s, e = row['consensus_abs_start'], row['consensus_abs_end']
    matching_segment = df[(df['abs_end'] >= s) & (df['abs_start'] <= e)]
    if len(matching_segment) > 1:
        matching_segment['overlap'] = matching_segment.apply(lambda row: get_overlap(row, s, e), axis=1)
        matching_segment = matching_segment[matching_segment['overlap'] > 1]
    if len(matching_segment) == 1:
        return np.unique(matching_segment.segment.values)[0]


def get_overlap(row, s, e):
    x = range(s, e + 1)
    y = range(row['abs_start'], row['abs_end'] + 1)
    return len(range(max(x[0], y[0]), min(x[-1], y[-1]) + 1)) - 1


def columns_cleanup(df):
    # keep chr in seg name:
    common = ['consensus_abs_start', 'consensus_abs_end', 'consensus_segment']
    if len(df['consensus_segment'].iloc[0].split('_')) == 2:
        df['consensus_segment'] = df['segment'].str.split('_').str[0] + '_' + df['consensus_segment']
    A_cols = [c for c in df.columns if (('_A' in c) and ('_AB' not in c))]
    B_cols = [c for c in df.columns if ('_B' in c)]
    AB_cols = [c for c in df.columns if ('_AB' in c)]
    clone_cols = A_cols + B_cols + AB_cols
    df = df[common + clone_cols].drop_duplicates()
    df_a = pd.melt(df[common + A_cols], id_vars=common, value_vars=A_cols, var_name='clone', value_name='cn')
    df_a['clone'] = df_a['clone'].str.split('_', expand=True)[0]
    df_a['allele'] = 'A'
    df_b = pd.melt(df[common + B_cols], id_vars=common, value_vars=B_cols, var_name='clone', value_name='cn')
    df_b['clone'] = df_b['clone'].str.split('_', expand=True)[0]
    df_b['allele'] = 'B'
    try:
        df_ab = pd.melt(df[common + AB_cols], id_vars=common, value_vars=AB_cols, var_name='clone', value_name='cn')
        df_ab['clone'] = df_ab['clone'].str.split('_', expand=True)[0]
        df_ab['allele'] = 'AB'
    except:
        df_ab = pd.DataFrame(columns=df_a.columns)
    df = pd.concat([df_a, df_b, df_ab])
    return df


def sort_alleles(df):
    A = df[df['allele'] == 'A'].rename(columns={'cn': 'cn_L'}).drop(columns='allele')
    B = df[df['allele'] == 'B'].rename(columns={'cn': 'cn_R'}).drop(columns='allele')
    AB = df[df['allele'] == 'AB']
    d = pd.merge(A, B)
    mismatched = d[d['cn_L'] < d['cn_R']]
    mismatched_L = list(mismatched.cn_L)
    mismatched_R = list(mismatched.cn_R)
    row_indexer = d['cn_L'] < d['cn_R']
    d.loc[row_indexer, 'cn_L'] = mismatched_R
    d.loc[row_indexer, 'cn_R'] = mismatched_L
    A = d.drop(columns='cn_R').rename(columns={'cn_L': 'cn'})
    A['allele'] = 'A'
    B = d.drop(columns='cn_L').rename(columns={'cn_R': 'cn'})
    B['allele'] = 'B'
    df = pd.concat([A, B, AB])
    return df


def convert_alpaca_to_wide_format(a):
    a = a[['tumour_id','clone', 'pred_CN_A', 'pred_CN_B', 'segment']].drop_duplicates()
    a.rename(columns={'pred_CN_A': 'A', 'pred_CN_B': 'B'}, inplace=True)
    a = a.pivot(index='segment', columns='clone', values=['A', 'B'])
    new_index = [f'{x[1]}_{x[0]}' for x in a.columns]
    a.columns = new_index
    a.reset_index(inplace=True)
    a['chr'] = 'chr' + a['segment'].str.split('_', expand=True)[0]
    a['start'] = a['segment'].str.split('_', expand=True)[1].astype(int)
    a['end'] = a['segment'].str.split('_', expand=True)[2].astype(int)
    return a


def matching_comparison_with_true(df,true):
    def matching_comparison(seg_df):
        segment = seg_df.consensus_segment.unique()[0]
        seg_df['seg_len'] = seg_df.consensus_abs_end-seg_df.consensus_abs_start
        true_seg_df = true[true.consensus_segment == segment]
        if len(true_seg_df)==0:
            return pd.DataFrame()
        true_states_A = set(true_seg_df.loc[true_seg_df.allele=='A',['clone','cn']].cn.values)
        true_states_B = set(true_seg_df.loc[true_seg_df.allele=='B',['clone','cn']].cn.values)
        found_states_A = set(seg_df.loc[seg_df.allele=='A',['clone','cn']].cn.values)
        found_states_B = set(seg_df.loc[seg_df.allele=='B',['clone','cn']].cn.values)
        # how many of the true states were found?
        found_A = len(true_states_A.intersection(found_states_A))/len(true_states_A)
        found_B = len(true_states_B.intersection(found_states_B))/len(true_states_B)
        seg_df['found'] = 0
        seg_df.loc[seg_df.allele=='A','found'] = found_A
        seg_df.loc[seg_df.allele=='B','found'] = found_B
        return seg_df
    scored_df = df.groupby(['consensus_segment','clone']).apply(lambda x: matching_comparison(x)).reset_index(drop=True)
    scored_df = scored_df[['consensus_segment','allele','seg_len','found']].drop_duplicates()
    return scored_df


parser = argparse.ArgumentParser()
parser.add_argument('--input_data_directory',type=str,required=True,help='Directory where input data is stored. Should contain subdirectories for each tumour')
parser.add_argument('--case_results_file', type=str, required=True, help='Results for one tumour')
parser.add_argument('--chr_table_file', type=str, required=True, help='')

args = parser.parse_args()
input_data_directory = args.input_data_directory
alpaca = pd.read_csv(args.case_results_file)
chr_table_file = args.chr_table_file

tumour_id = alpaca.tumour_id.unique()[0]
tumour_dir = f'{input_data_directory}/{tumour_id}'
true = add_abs_coordinates(pd.read_csv(f'{tumour_dir}/copynumbers.csv'), chr_table_file)
hatchet = add_abs_coordinates(pd.read_csv(f'{tumour_dir}/hatchet.csv'), chr_table_file)
hd = add_abs_coordinates(pd.read_csv(f'{tumour_dir}/clonehd.csv'), chr_table_file)
alpaca = add_abs_coordinates(convert_alpaca_to_wide_format(alpaca), chr_table_file)

dfs = find_consensus_segmentation([alpaca, true,hatchet,hd])
dfs = [columns_cleanup(df) for df in dfs]
dfs = [sort_alleles(df) for df in dfs]
alpaca, true,hatchet,hd = dfs
all_segments = [set(alpaca.consensus_segment.unique()),set(true.consensus_segment.unique()),set(hatchet.consensus_segment.unique()),set(hd.consensus_segment.unique())]
common_segments = set.intersection(*all_segments)
unique_segments = set().union(*(s - common_segments for s in all_segments))
# keep only common segments:
alpaca = alpaca[alpaca.consensus_segment.isin(common_segments)]
true = true[true.consensus_segment.isin(common_segments)]
hatchet = hatchet[hatchet.consensus_segment.isin(common_segments)]
hd = hd[hd.consensus_segment.isin(common_segments)]

hatchet_vs_true = matching_comparison_with_true(hatchet, true.copy())
alpaca_vs_true = matching_comparison_with_true(alpaca, true.copy())
hd_vs_true = matching_comparison_with_true(hd, true.copy())
hatchet_vs_true_weighted = np.average(a=hatchet_vs_true['found'],weights=hatchet_vs_true['seg_len'])
alpaca_vs_true_weighted = np.average(a=alpaca_vs_true['found'],weights=alpaca_vs_true['seg_len'])
hd_vs_true_weighted = np.average(a=hd_vs_true['found'],weights=hd_vs_true['seg_len'])

results = pd.DataFrame({'tumour_id':[tumour_id],
                        'alpaca':[alpaca_vs_true_weighted],    
                        'hatchet':[hatchet_vs_true_weighted],    
                        'hd':[hd_vs_true_weighted]})
results.to_csv(f'{tumour_id}_matched_comparison.csv',index=False)
