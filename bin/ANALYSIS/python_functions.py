import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

def get_tree_edges(tree_paths):
    all_edges = list()
    for path in tree_paths:
        if len(path) == 2:
            all_edges.append(tuple(path))
        else:
            for i in range(len(path) - 1):
                all_edges.append((path[i], path[i + 1]))
    unique_edges = set(all_edges)
    return unique_edges


def get_chr_table(chr_table_path):
    chr_table = pd.read_csv(chr_table_path)
    chr_table['cumsum'] = np.cumsum(chr_table['len'])
    chr_table['shift'] = [0] + list(chr_table['cumsum'][:-1])
    chr_table['ticks'] = chr_table['shift'] + chr_table['len'] / 2
    chr_table = chr_table[:-2]
    return chr_table

def flat_list(target_list):
    if isinstance(target_list[0], list):
        return [item for sublist in target_list for item in sublist]
    else:
        return target_list
    
def clean_output(patient_output):
    if ('abs_segment' in patient_output.columns ) and ('segment' not in patient_output.columns):
        patient_output = patient_output.T.drop_duplicates().T
        patient_output['segment'] = patient_output['abs_segment']
    if len(patient_output["segment"].iloc[0].split('_')) == 3:
        patient_output["chr"] = patient_output["segment"].str.split("_", expand=True)[0].str.replace('chr', '').apply(lambda x: clean_segment(x))
        patient_output['start'] = patient_output["segment"].str.split("_", expand=True)[1].apply(lambda x: clean_segment(x))
        patient_output['end'] = patient_output["segment"].str.split("_", expand=True)[2].apply(lambda x: clean_segment(x))
    else:  # this means that chr column is already unpacked and that segment is in format 123_345 (no chr)
        patient_output['start'] = patient_output["segment"].str.split("_", expand=True)[0].apply(lambda x: clean_segment(x))
        patient_output['end'] = patient_output["segment"].str.split("_", expand=True)[1].apply(lambda x: clean_segment(x))
    patient_output = patient_output.loc[~(patient_output["chr"].isna() | (patient_output["start"].isna()) | (patient_output["end"].isna()))]
    return patient_output


def clean_segment(seg):
    if isinstance(seg, str):
        if seg.isdigit():
            seg = int(seg)
        else:
            seg = np.nan
    return seg


def read_tree_json(json_path):
    with open(json_path, 'r') as f:
        tree = json.load(f)
    return tree


def compare_alpaca_to_true(alpaca, true, chr_table_file, sort_alleles=False):
    alpaca = convert_alpaca_to_wide_format(alpaca)
    alpaca = add_abs_coordinates(alpaca, chr_table_file)
    true = add_abs_coordinates(true, chr_table_file)
    
    dfs = find_consensus_segmentation([alpaca, true])
    dfs = [columns_cleanup(df) for df in dfs]
    if sort_alleles:
        dfs = [sort_alleles(df) for df in dfs]
    
    alpaca = dfs[0]
    true = dfs[1]
    
    alpaca = pd.merge(alpaca, true, on=['consensus_abs_start', 'consensus_abs_end', 'consensus_segment', 'clone', 'allele'], how='inner', suffixes=('_alpaca', '_true'))
    alpaca['len'] = alpaca['consensus_abs_end'] - alpaca['consensus_abs_start']
    alpaca['correct'] = alpaca['cn_alpaca'] == alpaca['cn_true']
    alpaca['abs_segment'] = alpaca['consensus_segment']
    alpaca = alpaca.merge(convert_absolute_to_relative(alpaca['consensus_segment'], chr_table_file), on='consensus_segment')
    return alpaca.drop_duplicates()

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
    segments_are_the_same = [set(list_of_dfs[0].segment)==set(list_of_dfs[i].segment) for i in range(1,len(list_of_dfs))]
    # chcek if list contains only true values:
    if False in segments_are_the_same:
        print('Finding consensus segmentation')
        output_dfs = []
        # breakpoint is defined as 'ends', excluding the last end:
        new_breakpoints = list(set(flat_list([(list(x['abs_start']) + list(x['abs_end'])) for x in list_of_dfs])))
        new_breakpoints = sorted(new_breakpoints)
        for i, d in enumerate(list_of_dfs):
            d = d.sort_values('abs_start', ascending=True)
            new_starts = new_breakpoints[:-1]  # [x+1 for x in new_ends]
            new_ends = new_breakpoints[1:]  # [b for b in new_breakpoints if b not in list(d['abs_end'])]
            # new_coordinates = sorted(list(d.abs_start) +  list(d.abs_end) + new_ends + new_starts)
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
    else:
        #print('segments are the same')
        output_dfs = []
        for i, d in enumerate(list_of_dfs):
            d = d.sort_values('abs_start', ascending=True)
            d['consensus_abs_start'] = d['abs_start']
            d['consensus_abs_end'] = d['abs_end']
            d['consensus_segment'] = d['chr'].str.replace('chr', '') + '_' + d['abs_start'].astype(str) + '_' + d['abs_end'].astype(str)
            output_dfs.append(d)
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
    a['AB'] = a['A'] + a['B']
    a = a.pivot(index='segment', columns='clone', values=['A', 'B', 'AB'])
    new_index = [f'{x[1]}_{x[0]}' for x in a.columns]
    a.columns = new_index
    a.reset_index(inplace=True)
    a['chr'] = 'chr' + a['segment'].str.split('_', expand=True)[0]
    a['start'] = a['segment'].str.split('_', expand=True)[1].astype(int)
    a['end'] = a['segment'].str.split('_', expand=True)[2].astype(int)
    return a


def convert_absolute_to_relative(segment_series,chr_table_path):
    chr_table = get_chr_table(chr_table_path)
    df = pd.DataFrame(segment_series)
    df['chr'] = 'chr'+df['consensus_segment'].str.split('_', expand=True)[0]
    df['start_abs'] = df['consensus_segment'].str.split('_', expand=True)[1].astype(int)
    df['end_end'] = df['consensus_segment'].str.split('_', expand=True)[2].astype(int)
    df = df.merge(chr_table, on='chr')
    df['start'] = df['start_abs'] - df['shift']
    df['end'] = df['end_end'] - df['shift']
    df['segment'] = df['chr'].str.replace('chr','') + '_' + df['start'].astype(str) + '_' + df['end'].astype(str)
    return df[['consensus_segment','segment']]


def add_abs_coordinates(df,chr_table_path):
    chr_table = get_chr_table(chr_table_path)
    df = df.merge(chr_table, on='chr')
    df['abs_start'] = df['start'] + df['shift']
    df['abs_end'] = df['end'] + df['shift']
    df.drop(columns=['shift','len','cumsum','ticks'], inplace=True)
    return df

def get_tree_edges(tree_paths):
    all_edges = list()
    for path in tree_paths:
        if len(path) == 2:
            all_edges.append(tuple(path))
        else:
            for i in range(len(path) - 1):
                all_edges.append((path[i], path[i + 1]))
    unique_edges = set(all_edges)
    return unique_edges

def get_unique_lists(list_of_lists):
    u = []
    for e in list_of_lists:
        if e not in u:
            u.append(e)
    return u

def get_ancestors(clone, tree):
    branch = [b for b in tree if clone in b][0]
    ancestors = branch[:branch.index(clone) + 1]
    return ancestors

def get_descendants(clone, tree):
    descendants = []
    for branch in tree:
        if clone in branch:
            descendants.append(branch[branch.index(clone) + 1:])
    return list(set(flat_list(descendants)))


def find_parent(clone, tree):
    for branch in tree:
        if clone in branch:
            clone_index = branch.index(clone)
            if clone_index != 0:
                return branch[clone_index - 1]
            else:
                return 'diploid'
            
def get_edge(tree: list):
    edges = []
    tree_edges = []
    for branch in tree:
        branch_edges = []
        for i, c in enumerate(branch[:-1]):
            edges.append((branch[i], branch[i + 1]))
            branch_edges.append((branch[i], branch[i + 1]))
        tree_edges.append(branch_edges)
    edges = list(set(edges))
    return edges, tree_edges

def remove_duplicates_preserve_order(input_list):
    seen = {}
    result = []
    for item in input_list:
        if item not in seen:
            seen[item] = True
            result.append(item)
    return result


def get_driver_mutations(mut_table,tumour_id,chr_table_file):
    mutations = pd.read_csv(mut_table)
    mutations = mutations[mutations.tumour_id == tumour_id]
    driver_mutations = mutations[(mutations.DriverMut == True)][['tumour_id','chr','start','Hugo_Symbol','PyCloneCluster_SC']]
    driver_mutations['chr'] = 'chr' + driver_mutations['chr'].astype(str)
    driver_mutations = driver_mutations.merge(get_chr_table(chr_table_file)[['chr','shift']])
    driver_mutations['abs_position'] = driver_mutations['start'] + driver_mutations['shift']
    driver_mutations['clone'] = 'clone' + driver_mutations['PyCloneCluster_SC'].astype(str)
    driver_mutations['gene'] = driver_mutations['Hugo_Symbol']
    return driver_mutations[['tumour_id','abs_position','clone','gene']]


def create_seaborn_colour_map(linSegCmap, number_of_colours):
    colors = [linSegCmap(i) for i in np.linspace(0, 1, number_of_colours)]
    sns_palette = sns.color_palette(colors)
    return sns_palette

def get_colourmap(max_val,testing=False,cn_type='total'):
    min_val=0
    positive_colors = np.array([1,0,0,1])
    negative_colors = np.array([0,0,1,1])
    white = np.array([1,1,1,1])
    if cn_type == 'total':
        colors = [
            (0.0, negative_colors),  # Blue at 0
            (1.0 / max_val, negative_colors),  # Blue at 1
            (2.0 / max_val, white),  # White at 2
            (1.0, positive_colors)   # Red at max_val
            ]
    else:
        colors = [
            (0.0, negative_colors),  # Blue at 0
            (1.0 / max_val, white),  # White at 2
            (1.0, positive_colors)   # Red at max_val
            ]
    custom_colormap = LinearSegmentedColormap.from_list('custom_map', colors)
    # testing:
    if testing:
        data = np.random.uniform(low=min_val, high=max_val, size=(10,10))
        data = np.array(range(min_val, max_val + 1)).reshape(1, -1)
        plt.imshow(data, cmap=custom_colormap)
        plt.colorbar()
        plt.show()
    return custom_colormap