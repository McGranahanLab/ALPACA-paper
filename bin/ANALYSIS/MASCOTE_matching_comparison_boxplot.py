#!/usr/bin/env python3

import scipy.stats as stats
import pandas as pd
import numpy as np
import argparse
import plotly.graph_objects as go


def get_hd_hatchet_alpaca_colours(palette):
    import json
    with open(palette, 'r') as f:
        palette = json.load(f)
    catgorical_palette = palette['categorical']
    colours = {}
    colours['hd'] = catgorical_palette['c12']
    colours['hatchet'] = catgorical_palette['c8']
    colours['alpaca'] = catgorical_palette['c2']
    return colours


def pval_to_asterisk(p_val):
    if p_val < 0.001:
        return '***'
    elif p_val < 0.01:
        return '**'
    elif p_val < 0.05:
        return '*'
    else:
        return 'ns'


def effect_size_rank_biserial_and_p_val(x, y):
    # Perform Wilcoxon test
    w, p = stats.wilcoxon(x, y)
    # Rank the data
    ranks = stats.rankdata(np.concatenate((x, y)))
    # Number of elements in each group
    n1 = len(x)
    n2 = len(y)
    # Sums of ranks in each group
    R1 = np.sum(ranks[:n1])
    R2 = np.sum(ranks[n1:])
    # Rank-biserial correlation
    rank_biserial = 1 - 2 * (min(R1, R2)) / (n1 * n2)
    # r effect size. Not using w/np.sqrt(n1) because for clone hd w statistic is zero (clone hd is worse than alpaca in every case)
    differences = x - y
    ranks = stats.rankdata(np.abs(differences))
    positive_rank_sum = np.sum(ranks[differences > 0])
    negative_rank_sum = np.sum(ranks[differences < 0])
    # Use the smaller of the two rank sums
    R_smaller = min(positive_rank_sum, negative_rank_sum)
    r = 1 - (R_smaller / (n1 * (n1 + 1) / 2))
    return r, p, rank_biserial


def make_matching_comparison_plot(cohort_results, palette_path='../../_assets/publication_palette.json'):
    colours = get_hd_hatchet_alpaca_colours(palette_path)
    effect_size_hatchet, p_value_hatchet, rank_biserial_hatchet = effect_size_rank_biserial_and_p_val(cohort_results.alpaca, cohort_results.hatchet)
    effect_size_clone_hd, p_value_clone_hd, rank_biserial_hd = effect_size_rank_biserial_and_p_val(cohort_results.alpaca, cohort_results.hd)
    results = pd.DataFrame(
        {'method': ['HATCHet', 'cloneHD', 'ALPACA'],
         'effect_size': [effect_size_hatchet, effect_size_clone_hd, ''],
         'p_value_vs_alpaca': [p_value_hatchet, p_value_clone_hd, ''],
         'rank_biserial': [rank_biserial_hatchet, rank_biserial_hd, ''],
         'mean': [np.mean(cohort_results.hatchet), np.mean(cohort_results.hd), np.mean(cohort_results.alpaca)],
         'median': [np.median(cohort_results.hatchet), np.median(cohort_results.hd), np.median(cohort_results.alpaca)], })
    scale = 0.2
    font_size = 45
    w = 600
    h = 1600
    fonts = 'DejaVu Sans'
    fig = go.Figure()
    fig.update_layout(title='', font=dict(family=fonts, size=font_size), width=w, height=h, showlegend=False, plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)',
                      annotations=[
                          {'text': f'{pval_to_asterisk(p_value_hatchet)}', 'x': 0.75, 'y': 1.02, 'xref': 'paper', 'yref': 'paper', 'align': 'left', 'showarrow': False, 'font': {'size': 95}},
                          {'text': f'{pval_to_asterisk(p_value_clone_hd)}', 'x': 0.5, 'y': 1.08, 'xref': 'paper', 'yref': 'paper', 'align': 'left', 'showarrow': False, 'font': {'size': 95}}, ], )
    fig.update_xaxes(showgrid=False)
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='black', zeroline=True, zerolinecolor='black', zerolinewidth=1, title='Weighted matched scores')
    jitter_value = 0
    point_size = 7
    line_width = 2
    alpaca_colour = colours['alpaca']
    hatchet_colour = colours['hatchet']
    clone_hd_colour = colours['hd']
    clone_hd_trace = go.Box(y=cohort_results.hd, orientation='v', boxpoints='all', jitter=jitter_value, pointpos=0, line=dict(color=clone_hd_colour, width=line_width), marker=dict(size=point_size),
                            name='cloneHD')
    hatchet_trace = go.Box(y=cohort_results.hatchet, orientation='v', boxpoints='all', jitter=jitter_value, pointpos=0, line=dict(color=hatchet_colour, width=line_width), marker=dict(size=point_size),
                           name='HATCHet')
    alpaca_trace = go.Box(y=cohort_results.alpaca, orientation='v', boxpoints='all', jitter=jitter_value, pointpos=0, line=dict(color=alpaca_colour, width=line_width), marker=dict(size=point_size),
                          name='ALPACA')
    fig.add_trace(clone_hd_trace)
    fig.add_trace(hatchet_trace)
    fig.add_trace(alpaca_trace)
    fig.update_layout(
        paper_bgcolor='white',
        plot_bgcolor='white'
    )
    return fig, w, h, scale, results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scores_matched_comparison_all_tumours', type=str, required=True, help='Results for all tumours')
    parser.add_argument('--palette_path', type=str, required=True, help='Colour palette path')
    args = parser.parse_args()
    scores_matched_comparison_all_tumours = args.scores_matched_comparison_all_tumours
    scores_matched_comparison_all_tumours = pd.read_csv(args.scores_matched_comparison_all_tumours)
    pallete_path = args.pallete_path
    fig, w, h, scale, results = make_matching_comparison_plot(scores_matched_comparison_all_tumours, pallete_path)
    fig.write_image(f'matching_comparison_boxplot.png', width=w, height=h, scale=scale)
    results.to_csv('matching_comparison_p_vals.csv', index=False)


if __name__ == "__main__":
    main()
