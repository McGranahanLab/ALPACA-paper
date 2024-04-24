#!/usr/bin/env python3
import pandas as pd
import argparse
import matplotlib.pyplot as plt

def make_histogram(fig, axes,scores,metric,ax_number,colors,title,font_scale,set_y_ticks,font_name):
    input = scores[metric]
    ax = axes[ax_number]
    if 'Alpaca' in metric:
        color = colors['alpaca']
    else:
        color = colors['simple']
    input.hist(ax=ax, color=color)
    mean_value = input.mean()
    ax.axvline(mean_value, color='k', linestyle='--')
    ax.text(mean_value, ax.get_ylim()[1] * 0.7, f'Mean:\n {mean_value:.2f}', fontsize=8 * font_scale)
    ax.set_title(f'{metric}', fontname=font_name, fontsize=12 * font_scale)
    ax.set_facecolor('white')
    if set_y_ticks:
        ax.set_yticks([0,500,1000,1500,2000])
    for tick in ax.get_xticklabels():
        tick.set_fontname(font_name)
        tick.set_fontsize(8 * font_scale)
    for tick in ax.get_yticklabels():
        tick.set_fontname(font_name)
        tick.set_fontsize(8 * font_scale)


def make_comparison_histogram(alpaca_scores,simple_scores,alpaca_colour = '#1F78B4',simple_colour = '#FB9A99',set_y_ticks = True,font_name = 'DejaVu Sans'):
    alpaca_scores = alpaca_scores[['tumour_id','clone','clone_precision','clone_recall']].copy()
    simple_scores = simple_scores[['tumour_id','clone','clone_precision','clone_recall']].copy()
    alpaca_scores.rename(columns={'clone_precision':'Alpaca Precision'}, inplace=True)
    alpaca_scores.rename(columns={'clone_recall':'Alpaca Recall'}, inplace=True)
    simple_scores.rename(columns={'clone_precision':'Simple Model Precision'}, inplace=True)
    simple_scores.rename(columns={'clone_recall':'Simple Model Recall'}, inplace=True)
    scores = pd.merge(alpaca_scores, simple_scores, on=['tumour_id','clone'])
    scores = scores[['tumour_id','clone','Alpaca Precision','Simple Model Precision','Alpaca Recall','Simple Model Recall']]
    fig, axes = plt.subplots(4, 1, figsize=(10, 15), sharex=True)
    metrics = ['Alpaca Precision','Simple Model Precision','Alpaca Recall','Simple Model Recall']
    font_scale = 3
    if len(scores) > 0:
        for j, metric in enumerate(metrics):
            make_histogram(fig, axes,scores,metric,ax_number=j,colors={'alpaca':alpaca_colour,'simple':simple_colour},title="",font_scale=font_scale,set_y_ticks=set_y_ticks,font_name=font_name)
    fig.subplots_adjust(left=0.2)
    fig.text(0.5, 0.04, 'Score per clone', ha='center', va='center', fontname=font_name, fontsize=12 * font_scale)
    fig.text(0.04, 0.5, 'Frequency', ha='center', va='center', rotation='vertical', fontname=font_name, fontsize=12 * font_scale)
    return fig
    
            
def main():        
    parser = argparse.ArgumentParser()
    parser.add_argument('--alpaca_scores', type=str, required=True, help='Alpaca results scord against known copy-numbers')
    parser.add_argument('--alpaca_colour', type=str, required=False, default='#1F78B4', help='Colour string for the histogram') 
    parser.add_argument('--alpaca_name', type=str, required=False, default='ALPACA')
    parser.add_argument('--simple_scores', type=str, required=True, help='Simple model results scord against known copy-numbers')
    parser.add_argument('--simple_colour', type=str, required=False, default='#FB9A99', help='Colour string for the histogram')
    parser.add_argument('--simple_name', type=str, required=False, default='Simple Model')
    parser.add_argument('--set_y_ticks', type=bool, required=False, default=False)
    print('exe')
    args = parser.parse_args()
    alpaca_scores = pd.read_csv(args.alpaca_scores)
    alpaca_colour = args.alpaca_colour
    alpaca_name = args.alpaca_name

    simple_scores = pd.read_csv(args.simple_scores)
    simple_colour = args.simple_colour
    simple_name = args.simple_name
    set_y_ticks = args.set_y_ticks
    font_name = 'DejaVu Sans'

    '''
    #TODO remove DEV
    alpaca_scores = pd.read_csv('../output/simulation/simulations_new_constr/cohort_outputs/subclonal_scores_summary_ALPACA_simulations_new_constr.csv')
    simple_scores = pd.read_csv('../output/simulation/simulations_new_constr/cohort_outputs/subclonal_scores_summary_ALPACA_simulations_new_constr.csv')
    alpaca_colour = '#1F78B4'
    simple_colour = '#FB9A99'
    set_y_ticks = True
    '''

    fig = make_comparison_histogram(alpaca_scores,simple_scores)
    fig.savefig(f'histograms_alpaca_simple_comparison.png', bbox_inches='tight')

if __name__ == "__main__":
    main()