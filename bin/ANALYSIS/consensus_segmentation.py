import pandas as pd
import numpy as np


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
    assert (len(set([x['consensus_segment'].nunique() for x in output_dfs])) == 1)
    return output_dfs


def find_original_matching_segment(row, df):
    s, e = row['consensus_abs_start'], row['consensus_abs_end']
    matching_segment = df[(df['abs_end'] >= s) & (df['abs_start'] <= e)]
    matching_segment = matching_segment[['chr','abs_start','abs_end','segment']].drop_duplicates()
    if len(matching_segment) > 1:
        matching_segment['overlap'] = matching_segment.apply(lambda matching_segment_row: get_overlap(matching_segment_row, s, e), axis=1)
        matching_segment = matching_segment[matching_segment['overlap'] > 1]
    if len(matching_segment) == 1:
        return np.unique(matching_segment.segment.values)[0]
    

def get_overlap(row, s, e):
    x = range(s, e + 1)
    y = range(row['abs_start'], row['abs_end'] + 1)
    return len(range(max(x[0], y[0]), min(x[-1], y[-1]) + 1)) - 1


def flat_list(target_list):
    if isinstance(target_list[0], list):
        return [item for sublist in target_list for item in sublist]
    else:
        return target_list


def get_chr_table(chr_table_path):
    chr_table = pd.read_csv(chr_table_path)
    chr_table['cumsum'] = np.cumsum(chr_table['len'])
    chr_table['shift'] = [0] + list(chr_table['cumsum'][:-1])
    chr_table['ticks'] = chr_table['shift'] + chr_table['len'] / 2
    chr_table = chr_table[:-2]
    return chr_table


def add_abs_coordinates(df,chr_table_path):
    chr_table = get_chr_table(chr_table_path)
    df = df.merge(chr_table, on='chr')
    df['abs_start'] = df['start'] + df['shift']
    df['abs_end'] = df['end'] + df['shift']
    df.drop(columns=['shift','len','cumsum','ticks'], inplace=True)
    return df


def convert_absolute_to_relative(df,chr_table_path):
    df.drop(columns=['start','end','segment'], inplace=True)
    chr_table = get_chr_table(chr_table_path)
    df = df.merge(chr_table, on='chr')
    df['start'] = df['consensus_abs_start'] - df['shift']
    df['end'] = df['consensus_abs_end'] - df['shift']
    df['segment'] = df['chr'].str.replace('chr','') + '_' + df['start'].astype(str) + '_' + df['end'].astype(str)
    return df