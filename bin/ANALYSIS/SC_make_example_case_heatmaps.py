#!/usr/bin/env python3
import argparse
from call_wgd import call_wgd
import seaborn as sns
import plotly.graph_objects as go
import pandas as pd
import networkx as nx
from plotly.subplots import make_subplots
from python_functions import get_chr_table, clean_output, flat_list, read_tree_json, compare_alpaca_to_true, get_unique_lists

def true_unify_format(t):
    if 'clone_cpn' in t.columns:
        t['segment'] = t['segment'].str.replace('chr', '')
        pivoted_df = t.pivot(index='segment', columns=['clone', 'allele'], values='clone_cpn')
        new_columns = [f'{c[0]}_{c[1]}_true' for c in pivoted_df.columns]
        pivoted_df.columns = new_columns
        pivoted_df = pivoted_df.reset_index()
        pivoted_df['POS'] = pivoted_df['segment'].str.split('_', expand=True)[1]
        pivoted_df['CHR'] = pivoted_df['segment'].str.split('_', expand=True)[0]
        return pivoted_df
    else:
        return t


def get_colour(cp_state, color_map):
    if cp_state >= 8:
        return f'rgb{tuple([c * 255 for c in color_map[-1]])}'
    else:
        colours = {
            0: (0, 0, 1), 1: (0.66, 0.66, 0.66), 2: color_map[0], 3: color_map[1], 4: color_map[2], 5: color_map[3], 6: color_map[4],
            7: color_map[5]
        }
        return f'rgb{tuple([c * 255 for c in colours[cp_state]])}'


def plot_heat_map_with_true_solution(patient_output, allele, fig, tree_graph_df, color_map, number_of_clones, chr_table_file, plot_comparison=False):
    clones = tree_graph_df.sort_values('y_loc', ascending=True).clone
    clone_number = len(clones)
    chromosome_table = get_chr_table(chr_table_file)
    patient_output = clean_output(patient_output)
    patient_output['predicted_cpn'] = patient_output[f'pred_CN_{allele}']
    patient_output = patient_output.merge(tree_graph_df)
    patient_output = patient_output.sort_values('y_loc')
    color_map = sns.color_palette("rocket", 7, as_cmap=False)
    color_map.reverse()
    for clone_index, clone_name in enumerate(clones):
        clone_df = patient_output[patient_output['clone'] == clone_name]
        i = int(tree_graph_df[tree_graph_df.clone == clone_name]['y_loc'].iloc[0])
        i = list(reversed(range(number_of_clones)))[i]
        for cp_state in clone_df.predicted_cpn.unique():
            clone_df_cp_state = clone_df[clone_df.predicted_cpn == cp_state]
            segments_predicted = [[tuple([row[1]['abs_start'], 0]), tuple([row[1]['abs_end'], 1])] for row in clone_df_cp_state.iterrows()]
            segments_predicted_unique = []
            for x in segments_predicted:
                if x not in segments_predicted_unique:
                    segments_predicted_unique.append(x)
            cpn_color = get_colour(cp_state, color_map)
            for rectangle in segments_predicted_unique:
                segment = clone_df_cp_state[clone_df_cp_state['abs_start'] == rectangle[0][0]].segment.unique()[0]
                fig.add_trace(go.Scatter(
                    showlegend=False,
                    x=[rectangle[0][0], rectangle[0][0], rectangle[1][0], rectangle[1][0]],
                    y=[rectangle[0][1], rectangle[1][1], rectangle[1][1], rectangle[0][1]],
                    fill='toself',
                    mode='lines',
                    fillcolor=cpn_color,
                    line_color=cpn_color,
                    name=f'clone: {rectangle[0][1]}, seg: {segment}'), row=i + 1, col=2)
        chromosomes = [[tuple([row[1]['cumsum'], 0]), tuple([row[1]['cumsum'], clone_number])] for row in chromosome_table.iterrows()]
        for chromosome_line in chromosomes:
            fig.add_trace(go.Scatter(
                x=[chromosome_line[0][0], chromosome_line[0][0]],
                y=[0, 1],
                mode='lines', line=dict(color='black', width=1, dash='dash'), showlegend=False), row=i + 1, col=2)
    return fig


def plot_error_track(patient_output, allele, fig, tree_graph_df, number_of_clones, plot_black_markers=True, clones_to_skip=[]):
    if len(patient_output) == 0:
        return fig
    clones = tree_graph_df.sort_values('y_loc', ascending=True).clone
    patient_output = clean_output(patient_output)
    patient_output['predicted_cpn'] = patient_output[f'pred_CN_{allele}']
    patient_output = patient_output.merge(tree_graph_df)
    patient_output = patient_output.sort_values('y_loc')
    color_map = sns.color_palette("rocket", 7, as_cmap=False)
    color_map.reverse()
    
    for clone_index, clone_name in enumerate(clones):
        if bool(clones_to_skip and (clone_name in clones_to_skip)):
            continue
        clone_df = patient_output[patient_output['clone'] == clone_name]
        i = int(tree_graph_df[tree_graph_df.clone == clone_name]['y_loc'].iloc[0])
        i = list(reversed(range(number_of_clones)))[i]
        if len(clone_df) > 0:
            for cp_state in clone_df.predicted_cpn.unique():
                clone_df_cp_state = clone_df[clone_df.predicted_cpn == cp_state]
                segments_predicted = [[tuple([row[1]['abs_start'], 0]), tuple([row[1]['abs_end'], 0.5])] for row in clone_df_cp_state.iterrows()]
                segments_predicted_unique = []
                for x in segments_predicted:
                    if x not in segments_predicted_unique:
                        segments_predicted_unique.append(x)
                cpn_color = get_colour(cp_state, color_map)
                for rectangle in segments_predicted_unique:
                    segment = clone_df_cp_state[clone_df_cp_state['abs_start'] == rectangle[0][0]].segment.unique()[0]
                    fig.add_trace(go.Scatter(
                        showlegend=False,
                        x=[rectangle[0][0], rectangle[0][0], rectangle[1][0], rectangle[1][0]],
                        y=[rectangle[0][1], rectangle[1][1], rectangle[1][1], rectangle[0][1]],
                        fill='toself',
                        mode='lines',
                        fillcolor=cpn_color,
                        line_color=cpn_color,
                        name=f'clone: {rectangle[0][1]}, seg: {segment}'), row=i + 1, col=2)
            if plot_black_markers:
                segments_error_x = flat_list([[v[1]['abs_start'], v[1]['abs_end'], None] for v in clone_df.iterrows()])
                segments_error_y = flat_list([[0.5, 0.5, None] for v in clone_df.iterrows()])
                fig.add_trace(go.Scatter(
                    showlegend=False,
                    hoverinfo="skip",
                    x=segments_error_x,
                    y=segments_error_y,
                    mode='lines',
                    fillcolor='black',
                    line=dict(
                        color='black',
                        width=1)),
                    row=i + 1, col=2)
    return fig


def coerce_to_alpaca_format(df, allele):
    df = df[df.allele == allele]
    df[f'pred_CN_{allele}'] = df['cn_true']
    df['abs_start'] = df['consensus_abs_start']
    df['abs_end'] = df['consensus_abs_end']
    df.rename(columns={'consensus_seg_name': 'abs_segment'}, inplace=True)
    df = df[[f'pred_CN_{allele}', 'abs_start', 'abs_end', 'abs_segment', 'clone']].drop_duplicates()
    return df


class plot_heatmap_with_tree_compare_with_true_solution_publication:
    def __init__(self, alpaca_output, input_data_directory, chr_table_file, wgd_clones, max_cpn_cap=99, allele='A', true_solution_df='', plot_comparison=False, sort_alleles=False,
                 figure_attributes=None, clones_to_skip=[]):
        self.max_cpn_cap = max_cpn_cap
        color_map = sns.color_palette("rocket", 7, as_cmap=False)
        color_map.reverse()
        # default settings:
        
        if figure_attributes is None:
            figure_attributes = {'font_size': 29,
                                 'wgd_annotation_font_size': 25,
                                 'pred_true_font_size': 15,
                                 'chromosome_labels_font_size': 15,
                                 'region_number_font_size': 13,
                                 'colorbar_font_size': 15,
                                 'annotation_font_size': 22,
                                 'width_mm': 620,
                                 'height_mm': 620 * 0.62,
                                 'legend_y_offset': -0.2,
                                 'cbar_y_offset': 0,
                                 'bottom_margin': 45,
                                 'clone_prop_annotation_legend_x': 0.85,
                                 'clone_prop_annotation_legend_y': -0.34,
                                 'seg_solved_legend_x': 1,
                                 'seg_solved_legend_y': -0.04,
                                 'incorrectly_legend_x': 0.99,
                                 'incorrectly_legend_y': -0.06,
                                 'correctly_legend_x': 1,
                                 'correctly_legend_y': -0.1, }
        for key, value in figure_attributes.items():
            globals()[key] = value
        DPI = 96
        h = int(height_mm * DPI / 25.4)
        w = int(width_mm * DPI / 25.4)
        # w=1600
        # h=max(1000, 75 * self.number_of_clones)
        pie_colors = ['#626CAD', '#DD5D46']
        json_path = f'{input_data_directory}/tree_paths.json'
        tree = read_tree_json(json_path)
        tree = [[c for c in b if c in ['trunk', 'JII', 'JIII', 'JIV']] for b in tree]
        seeding_clones = []
        mrca = tree[0][0]
        # prepare input:
        self.alpaca_output = alpaca_output.copy()
        self.alpaca_output = self.alpaca_output[['tumour_id', 'clone', 'pred_CN_A', 'pred_CN_B', 'segment']].drop_duplicates()
        self.alpaca_output['chr'] = 'chr' + self.alpaca_output.segment.str.split('_', expand=True)[0]
        self.alpaca_output['Start'] = self.alpaca_output.segment.str.split('_', expand=True)[1].astype(int)
        self.alpaca_output['End'] = self.alpaca_output.segment.str.split('_', expand=True)[2].astype(int)
        self.true_solution_df = true_solution_df
        cp_table = pd.read_csv(f'{input_data_directory}/cp_table.csv', index_col='clone')
        cp_table = cp_table.loc[['trunk', 'JII', 'JIII', 'JIV']]
        # modify segment positions to absolute:
        self.alpaca_output = self.alpaca_output.merge(get_chr_table(chr_table_file))
        self.alpaca_output['abs_start'] = self.alpaca_output['Start'] + self.alpaca_output['shift']
        self.alpaca_output['abs_end'] = self.alpaca_output['End'] + self.alpaca_output['shift']
        self.alpaca_output = self.alpaca_output.sort_values(['abs_start'], ascending=False)
        self.alpaca_output.loc[self.alpaca_output['pred_CN_A'] > self.max_cpn_cap, 'pred_CN_A'] = self.max_cpn_cap
        self.alpaca_output.loc[self.alpaca_output['pred_CN_B'] > self.max_cpn_cap, 'pred_CN_B'] = self.max_cpn_cap
        x_min = self.alpaca_output['abs_start'].min()
        x_max = self.alpaca_output['abs_start'].max()
        self.number_of_clones = len(self.alpaca_output.clone.unique())
        max_levels = max([len(b) for b in tree])
        tree_with_levels = [dict(zip(b, range(0, len(b)))) for b in tree]
        tree_with_levels = pd.concat([pd.DataFrame(tree_with_levels[x], index=[0]).transpose() for x, _ in enumerate(tree_with_levels)]).reset_index().drop_duplicates().rename(
            columns={'index': 'clone', 0: 'level'})
        # make empty y_loc df
        clone_y_location = dict(zip(self.alpaca_output.clone.unique(), range(0, self.number_of_clones)))
        clone_y_location = pd.DataFrame(clone_y_location, index=[0]).transpose().reset_index().rename(columns={'index': 'clone', 0: 'y_loc'})
        clone_y_location['y_loc'] = 100
        # find sections:
        # section is a part of a path that requires its own horizontal space on the final graph
        ori_tree = tree.copy()
        sections = [[tree[0][0]]]
        while len(flat_list(ori_tree)) > 1:
            for i, branch in enumerate(ori_tree):
                if (branch != []) and (set(flat_list(sections)) != set(flat_list(tree))):
                    branching_clones = list(pd.Series(flat_list(ori_tree)).value_counts()[pd.Series(flat_list(ori_tree)).value_counts() > 1].index)
                    if branching_clones == []:
                        branching_clones = [mrca]
                    section_start = max([branch.index(x) for x in branching_clones if x in branch]) + 1
                    section = branch[section_start:]
                    ori_tree[ori_tree.index(branch)] = branch[:section_start]
                    ori_tree = get_unique_lists(ori_tree)
                    if section != []:
                        sections.append(section)
        # order sections according to proximity
        # start with on of the longest paths
        section_termini = [x[-1] for x in sections]
        initial_node = [s for s in sections if len(s) == max([len(x) for x in sections])][0][-1]
        # simplify tree graph:
        simple_tree = [[clone for clone in branch if clone in section_termini] for branch in tree]
        edges =  (simple_tree)
        G = nx.Graph()
        for edge in edges:
            G.add_edge(edge[0], edge[1], weight=0)
        nodes = [initial_node] + [x for x in list(G.nodes) if x is not initial_node]
        distance_to_nodes = dict(nx.all_pairs_shortest_path_length(G))
        processed_nodes = [[s for s in sections if s[-1] == initial_node][0][-1]]
        while len(nodes) > 0:
            node = processed_nodes[-1]
            neighbours = pd.DataFrame(distance_to_nodes[node], index=['val']).transpose()
            neighbours.drop(inplace=True, index=processed_nodes)
            try:
                neighbours.drop(inplace=True, index=node)
            except KeyError:
                pass
            neighbours = neighbours.reset_index().rename(columns={'index': 'clone'}).merge(tree_with_levels, how='left', on='clone')
            neighbours['level'] = neighbours['level'].astype(int)
            if len(neighbours) > 0:
                # if there are to neighbours equally close, choose the one which has higher level (i.e. is deeper in the tree)
                # to do that, multiply proximity by negative level:
                closest_candidates = neighbours[neighbours.val == min(neighbours.val)]
                closest_neighbours = closest_candidates[closest_candidates.level == max(closest_candidates.level)]
                closest_neighbour = closest_neighbours.clone.values[0]
                if closest_neighbour not in processed_nodes:
                    processed_nodes.append(closest_neighbour)
            nodes = [n for n in nodes if n != node]
        sorted_sections_raw = []
        for node in processed_nodes:
            sorted_sections_raw.append([s for s in sections if s[-1] == node][0])
        # join single-clone sections below MRCA with their descendants:
        sorted_sections = []
        below_MRCA = True
        skip_element = False
        for ss in sorted_sections_raw:
            if skip_element:
                skip_element = False
                continue
            if mrca in ss:
                below_MRCA = False
            if below_MRCA:
                if len(ss) == 1:
                    joined = ss + sorted_sections_raw[sorted_sections_raw.index(ss) + 1]
                    sorted_sections.append(joined)
                    skip_element = True
                else:
                    sorted_sections.append(ss)
            else:
                sorted_sections.append(ss)
        # assign y location on the plot to sorted sections:
        available_y_locs = range(0, self.number_of_clones)
        below_MRCA = True
        for section in sorted_sections:
            if mrca in section:
                below_MRCA = False
            locs_for_this_section = list(available_y_locs[:len(section)])
            if below_MRCA:
                locs_for_this_section = list(reversed(locs_for_this_section))
            locs_for_this_section_dict = dict(zip(section, locs_for_this_section))
            available_y_locs = available_y_locs[len(section):]
            for n in locs_for_this_section_dict.keys():
                clone_y_location.loc[clone_y_location.clone == n, 'y_loc'] = locs_for_this_section_dict[n]
        tree_graph_df = pd.merge(tree_with_levels, clone_y_location)
        clone_prop_title = 'Clone proportions'
        s = [[{"type": "xy", "rowspan": self.number_of_clones}, {"type": "xy"}, {"type": "xy"}, {"type": "pie"}]]
        for c in range(self.number_of_clones - 1):
            s.append([None, {"type": "xy"}, {"type": "xy"}, {"type": "pie"}])
        self.fig = make_subplots(
            rows=self.number_of_clones, cols=4,
            column_widths=[0.1, 0.75, 0.1, 0.05],
            specs=s, horizontal_spacing=0.01, vertical_spacing=0.01,
            subplot_titles=('', '', '', ''))
        for clone_pos in tree_graph_df.y_loc:
            hline = go.Scatter(
                showlegend=False,
                x=[-0.3, max_levels],
                y=[clone_pos + 0.5, clone_pos + 0.5],
                mode='lines',
                line=dict(color='Green', dash='dot'))
            self.fig.append_trace(hline, row=1, col=1)
        # *** plot tree ***
        for branch in tree:
            branch_df = tree_graph_df[tree_graph_df.clone.isin(branch)]
            self.fig.append_trace(go.Scatter(
                showlegend=False,
                name='tree',
                x=branch_df['level'],
                y=branch_df['y_loc'],
                mode='lines+markers',
                marker=dict(
                    symbol='circle',
                    color='purple',
                    size=10,
                    line=dict(
                        color='purple',
                        width=2)),
                text=branch_df['clone']),
                row=1, col=1
            )
        if wgd_clones != []:
            gd_df = tree_graph_df[tree_graph_df.clone.isin(wgd_clones)]
            self.fig.append_trace(go.Scatter(
                showlegend=False,
                name='gd clone',
                x=gd_df['level'],
                y=gd_df['y_loc'],
                mode='markers+text',
                marker=dict(
                    symbol='circle-open',
                    color='rgb(255, 37, 62)',
                    size=16,
                    line=dict(
                        color='rgb(255, 37, 62)',
                        width=3)),
                text=['WGD'] * len(wgd_clones),
                textfont=dict(size=wgd_annotation_font_size),
                textposition="middle right"
            ),
                row=1, col=1
            )
        if seeding_clones != []:
            seeding_clones_df = tree_graph_df[tree_graph_df.clone.isin(seeding_clones)]
            self.fig.append_trace(go.Scatter(
                showlegend=False,
                name='seeding clone',
                x=seeding_clones_df['level'],
                y=seeding_clones_df['y_loc'],
                mode='markers+text',
                marker=dict(
                    symbol='circle-open',
                    color='steelblue',
                    size=19,
                    line=dict(
                        color='steelblue',
                        width=3)),
                text=['Seeding'],
                textfont=dict(size=wgd_annotation_font_size),
                textposition="bottom right"
            ),
                row=1, col=1
            )
        self.fig.update_yaxes(
            showgrid=False,
            tickmode='array',
            tickvals=list(tree_graph_df.sort_values('y_loc').y_loc),
            ticktext=[x.replace('clone', 'Clone ') for x in list(tree_graph_df.sort_values('y_loc').clone)],
            range=[-0.3, self.number_of_clones - 0.7],
            showticklabels=True,
            zeroline=False, row=1, col=1)
        self.fig.update_xaxes(
            showgrid=False, zeroline=False, row=1, col=1
        )
        # *** plot clones ***
        self.df = self.alpaca_output
        if self.true_solution_df is not None:
            self.alpaca_with_true = compare_alpaca_to_true(alpaca=self.alpaca_output, true=self.true_solution_df, chr_table_file=chr_table_file, sort_alleles=sort_alleles)
            self.df_true = coerce_to_alpaca_format(self.alpaca_with_true, allele)
            self.df_errors = coerce_to_alpaca_format(self.alpaca_with_true[self.alpaca_with_true.correct == False], allele)  # this has no chr column and '123_456' seg name
        else:
            self.df = self.alpaca_output
        y_limit = self.df[f'pred_CN_{allele}'].max()
        shapes = []
        chr_len = get_chr_table(chr_table_file)
        ######
        self.fig = plot_heat_map_with_true_solution(self.df.copy(), allele, self.fig, tree_graph_df, color_map, self.number_of_clones, chr_table_file, plot_comparison=False)
        if plot_comparison:
            # self.fig = plotHeatMap(self.df_true.copy(), allele, self.fig, tree_graph_df, color_map,self.number_of_clones,plot_comparison=plot_comparison) # plots with plot_comparison = true will be plotted in the lower half of each clone row
            # add error highlight
            self.fig = plot_error_track(self.df_errors, allele, self.fig, tree_graph_df, self.number_of_clones, plot_black_markers=True, clones_to_skip=clones_to_skip)
        # rename columns to keep just Rx instead of the full name
        # cp_table=cp_table.rename(columns={x:x.split('.R')[1] if 'LTX' in x else x for x in cp_table.columns})
        # cp_table.set_index('clone',inplace=True)
        # cp_table = cp_table[['1','2','3','4','5','6','7','8','9','10','11','12','13']]
        for clone in self.df.clone.unique():
            i = int(tree_graph_df[tree_graph_df.clone == clone]['y_loc'].iloc[0])
            i = list(reversed(range(self.number_of_clones)))[i]
            self.fig.update_xaxes(range=[x_min, x_max], showgrid=False, row=i + 1, col=2)
            # add dotted line dividing predicted from true:
            if clone not in clones_to_skip:
                self.fig.add_trace(go.Scatter(
                    showlegend=False,
                    hoverinfo="skip",
                    x=[0, self.alpaca_output['abs_start'].max()],
                    y=[0.5, 0.5],
                    mode='lines',
                    fillcolor='black',
                    line=dict(
                        dash='dot',
                        color='black',
                        width=0.5)),
                    row=i + 1, col=2)
            # add proportions in regions:
            clone_cp = cp_table.loc[[clone]]
            # clone_cp = clone_cp.transpose().reset_index()
            # clone_cp = clone_cp.loc[1:]
            # clone_cp.columns = ['regions', 'proportions']
            # cp_table.set_index('clones', inplace=True)
            showscale = i == len(self.df.clone.unique()) - 1
            clone_proportion_heatmap = go.Heatmap(
                z=clone_cp.values,
                x=clone_cp.columns,
                y=clone_cp.index,
                colorscale='Blues',
                showscale=showscale,
                colorbar=dict(
                    tickfont=dict(size=colorbar_font_size),
                    orientation='h',
                    x=0.89,
                    y=cbar_y_offset,
                    yanchor='top',
                    len=0.1,
                    thickness=20),
                hoverinfo='z',
                zauto=False,
                zmin=0,
                zmax=1,
            )
            self.fig.add_trace(clone_proportion_heatmap, row=i + 1, col=3)
            clone_comparison = self.alpaca_with_true[(self.alpaca_with_true.clone == clone) & (self.alpaca_with_true.allele == allele)].drop_duplicates()
            clone_correct = clone_comparison.correct.sum()
            clone_incorrect = clone_comparison.shape[0] - clone_correct
            labels = ['Correct', 'Incorrect']
            values = [clone_correct, clone_incorrect]
            pie_chart = go.Pie(labels=labels, values=values, textinfo='none', showlegend=False, marker=dict(colors=pie_colors))
            self.fig.add_trace(pie_chart, row=i + 1, col=4)
            if plot_comparison:
                if clone not in clones_to_skip:
                    tick_annotation_labels = ['true (t) ', 'predicted (p) '] if i == 0 else ['t ', 'p ']
                else:
                    tick_annotation_labels = ['', '']
            self.fig.update_yaxes(tickmode='array',
                                  tickvals=[0.25, 0.75],
                                  ticktext=tick_annotation_labels,
                                  tickfont_size=pred_true_font_size,
                                  showticklabels=True, row=i + 1, col=2)
            self.fig.update_yaxes(showticklabels=False, row=i + 1, col=3)
            if i != len(self.df.clone.unique()) - 1:
                self.fig.update_xaxes(showticklabels=False, row=i + 1, col=1)
                self.fig.update_xaxes(showticklabels=False, row=i + 1, col=2)
                self.fig.update_xaxes(showticklabels=False, row=i + 1, col=3)
            else:
                self.fig.update_xaxes(
                    tickmode='array',
                    tickfont=dict(size=chromosome_labels_font_size),
                    tickvals=chr_len['cumsum'] - (chr_len['len'] / 2),
                    ticktext=[str(x) for x in list(range(1, 23))]
                    , showticklabels=True, row=i + 1, col=2)
                # region number annotations:
                self.fig.update_xaxes(
                    tickmode='array',
                    tickangle=0,
                    tickfont=dict(size=region_number_font_size),
                    ticktext=clone_cp.columns
                    , showticklabels=True, row=i + 1, col=3)
        self.tree_graph_df = tree_graph_df
        # subtitle font size:
        self.fig.update_annotations(font_size=font_size)
        self.fig.update_layout(
            # title=f'{tumour_id}',
            plot_bgcolor='rgba(255,255,255,0)',
            autosize=False,
            width=w,
            height=h,
            legend_tracegroupgap=50,
            legend=dict(
                orientation='h',
                yanchor="top",
                y=-0.025,
                xanchor="left",
                x=0,
                font=dict(
                    size=font_size)
            ))
        # build legend:
        legend_row = 1
        self.fig.add_trace(go.Scatter(
            legendgroup='0',
            x=[None],
            y=[None],
            mode="markers",
            name="Whole Genome Doubling",
            marker=dict(size=100, color='rgb(255, 37, 62)', symbol='circle-open', line=dict(color='rgb(255, 37, 62)', width=50)),
        ), row=legend_row, col=1)
        # dummy legend to add annotation:
        self.fig.add_trace(go.Scatter(
            legendgroup='1',
            x=[None],
            y=[None],
            mode="markers",
            line=dict(
                color='white',
                width=3),
            name="     Copy number states:",
            marker=dict(size=10, color='white', symbol='circle-open', line=dict(color='white', width=3)),
        ), row=legend_row, col=1)
        legend = {
            '0': f'rgb{(0, 0, 255)}',
            '1': f'rgb{(168, 168, 168)}',
            '2': f'rgb{tuple([c * 255 for c in color_map[0]])}',
            '3': f'rgb{tuple([c * 255 for c in color_map[1]])}',
            '4': f'rgb{tuple([c * 255 for c in color_map[2]])}',
            '5': f'rgb{tuple([c * 255 for c in color_map[3]])}',
            '6': f'rgb{tuple([c * 255 for c in color_map[4]])}',
            '7': f'rgb{tuple([c * 255 for c in color_map[5]])}',
            '8+': f'rgb{tuple([c * 255 for c in color_map[-1]])}'
        }
        l_group = ['2', '3',
                   '4', '5',
                   '6', '7',
                   '8', '9', '10']
        for j, c in enumerate(legend.keys()):
            self.fig.add_trace(go.Scatter(
                legendgroup=l_group[j],
                x=[None],
                y=[None],
                mode='markers',
                name=f'{c}',
                marker=dict(
                    color=legend[c],
                    size=10, ),
                showlegend=True
            ), row=legend_row, col=1)
        annotations = [
            {'text': f'Clone proportions<br>in samples:', 'x': clone_prop_annotation_legend_x, 'y': clone_prop_annotation_legend_y, 'xref': 'paper', 'yref': 'paper', 'align': 'left',
             'showarrow': False, 'font': {'size': annotation_font_size}},
            {'text': f'Seg. solved:', 'x': seg_solved_legend_x, 'y': seg_solved_legend_y, 'xref': 'paper', 'yref': 'paper', 'align': 'left', 'showarrow': False,
             'font': {'size': annotation_font_size}},
            {'text': f'incorrectly', 'x': incorrectly_legend_x, 'y': incorrectly_legend_y, 'xref': 'paper', 'yref': 'paper', 'align': 'left', 'showarrow': False,
             'font': {'size': annotation_font_size, 'color': '#DD5D46'}},
            {'text': f'correctly', 'x': correctly_legend_x, 'y': correctly_legend_y, 'xref': 'paper', 'yref': 'paper', 'align': 'left', 'showarrow': False,
             'font': {'size': annotation_font_size, 'color': '#616CAD'}},
        ]
        '''if allele=='A':
            annotations.append({'text': f'16p:<br>predicted cn = 1<br>true cn = 2', 'x': 0.78, 'y': -0.08, 'xref': 'paper', 'yref': 'paper', 'align': 'left', 'showarrow': False,'font': {'size': annotation_font_size}})'''
        self.fig.update_layout(margin=dict(t=30, b=bottom_margin, l=10, r=20),
                               legend=dict(font=dict(size=font_size), y=legend_y_offset, yanchor='top'),
                               font=dict(family='Helvetica', size=font_size),
                               annotations=annotations
                               )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_data_directory", type=str, required=True, help='')
    parser.add_argument("--cohort_results_file", type=str, required=True, help='Results for all tumours in cohort')
    parser.add_argument("--selected_cases", type=str, required=True, help='')
    parser.add_argument("--chr_table_file", type=str, required=True, help='')
    parser.add_argument("--wgd_calls", type=str, required=True, help='')
    
    args = parser.parse_args()
    input_data_directory = args.input_data_directory
    cohort_results_file = args.cohort_results_file
    cohort_results = pd.read_csv(cohort_results_file)
    chr_table_file = args.chr_table_file
    wgd_calls = args.wgd_calls
    wgd_calls = pd.read_csv(wgd_calls)
    selected_tumours = args.selected_tumours
    selected_tumours = [x.strip() for x in selected_tumours[1:-1].split(',')]
    
    for tumour_id in selected_tumours:
        tumour_input_directory = f'{input_data_directory}/{tumour_id}'
        true = true_unify_format(pd.read_csv(f'{tumour_input_directory}/copynumbers.csv'))
        tumour_wgd_calls = wgd_calls[wgd_calls.tumour_id == tumour_id]
        wgd_clones = list(wgd_calls[wgd_calls.GD > 0].clones)
        alpaca_output = cohort_results[cohort_results.tumour_id == tumour_id]
        alpaca_output = alpaca_output[alpaca_output.clone.isin(['trunk', 'JII', 'JIII', 'JIV'])]
        for allele in ['A', 'B']:
            heatmap = plot_heatmap_with_tree_compare_with_true_solution_publication(alpaca_output=alpaca_output, input_data_directory=tumour_input_directory, chr_table_file=chr_table_file,
                                                                                    wgd_clones=wgd_clones, max_cpn_cap=8, allele=allele, true_solution_df=true, plot_comparison=True,
                                                                                    sort_alleles=False)
            heatmap.fig.write_image(f'{tumour_id}_{allele}_example_heatmap.pdf')


if __name__ == "__main__":
    main()
