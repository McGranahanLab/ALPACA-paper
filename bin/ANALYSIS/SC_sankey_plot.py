#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import itertools
import plotly.graph_objects as go
import argparse

'''
For publication, run the notebook SC_sankey_plot.ipynb and rearange the nodes in the interactive output
'''

def getRBrewerPairedPalette(alpha):
    p = ["#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"]
    rgb = []
    for h in p:
        rgb_colour = tuple(int(h.lstrip('#')[i:i + 2], 16) for i in (0, 2, 4))
        rgb_colour_with_alpha = f'rgba{(rgb_colour + (alpha,))}'
        rgb.append(rgb_colour_with_alpha)
    return rgb


def plotFlowDiagramSC(df, plotTitle,hide_labels):
    c_map_full = [  'rgba(31, 120, 180, 0.5)',
                    'rgba(168, 168, 168, 0.5)',
                    'rgba(178, 223, 138, 0.5)',
                    'rgba(251, 154, 153, 0.5)',
                    'rgba(227, 26, 28, 0.5)',
                    'rgba(227, 227, 227, 1)']
    unique_cn_true_states = np.sort(df.cn_true.unique())
    unique_cn_predicted_states = np.sort(df.cn_alpaca.unique())
    true_labels = [f"True cn {x}" for x in unique_cn_true_states]
    predicted_labels = [f"Predicted cn {x}" for x in unique_cn_predicted_states]
    all_labels = true_labels + predicted_labels
    if hide_labels==True:
        all_labels = ['' for x in all_labels]
    flow_pairs = list(itertools.product(unique_cn_true_states, unique_cn_predicted_states))
    fp = []
    fv = []
    sources = []
    targets = []
    for flow_pair in flow_pairs:
        s = flow_pair[0]
        t = flow_pair[1]
        sources.append(s)
        targets.append(t + len(unique_cn_true_states) + 1)
        flow_df = df[(df.cn_true == s) & (df.cn_alpaca == t)]
        flow_value = len(flow_df)
        fp.append(flow_pair)
        fv.append(flow_value)
    state_to_indx = dict(zip(
        sorted(list(set(sources))) +
        sorted(list(set(targets))),
        list(range(0, len(all_labels)))))
    sources = [state_to_indx[x] for x in sources]
    targets = [state_to_indx[x] for x in targets]
    c_map_plot = [c_map_full[x] for x in sources]
    fig = go.Figure(data=[go.Sankey(
        arrangement="snap",
        node=dict(
            pad=20,
            thickness=60,
            line=dict(color="black", width=1),
            label = all_labels,
            x = [0 for x in true_labels] + [1 for x in predicted_labels],
            y = [1/len(true_labels) * x for x in range(len(true_labels))]+[1/len(predicted_labels) * x for x in range(len(predicted_labels))],
            #y = [0,0.2,0.4,0.6,0.8,0,0.16,0.33,0.5,0.66],
            color="white"
        ),
        link=dict(
            source=sources,
            target=targets,
            value=fv,
            color=c_map_plot
        ))])
    fig.update_layout(title_text=plotTitle, font_size=28, width=600, height=1000, font_family="helvetica")
    return fig


def findSegmentsFullyCorrect(df):
    if sum(df.correct==True)/df.shape[0] == 1:
        return True
    else:
        return False
 

def findClonesFullyCorrect(df):
    if sum(df.correct==True)/df.shape[0] == 1:
        return True
    else:
        return False

 
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scores_all_tumours', type=str, required=True, help='Alpaca results scord against known copy-numbers')
    args = parser.parse_args()
    scores_all_tumours = args.scores_all_tumours
    
    scores_all_tumours = pd.read_csv(scores_all_tumours)
    scores_all_tumours = scores_all_tumours[['clone','allele','segment','tumour_id','cn_true','cn_alpaca','correct']].drop_duplicates()
    scores_all_tumours.cn_alpaca.clip(upper=5,inplace=True)

    c_map_full = getRBrewerPairedPalette(0.5)
    fig = plotFlowDiagramSC(scores_all_tumours,'',hide_labels=True)
    fig.write_image(f'sankey_plot_publication_no_labels.pdf')
    fig = plotFlowDiagramSC(scores_all_tumours,'',hide_labels=False)
    fig.write_image(f'sankey_plot_publication_with_labels.pdf')

    # calculate metrics:

    fully_correct_segments = scores_all_tumours.groupby(['segment','tumour_id']).apply(findSegmentsFullyCorrect).reset_index()
    fully_correct_segments = f'{fully_correct_segments[0].sum() / fully_correct_segments.shape[0]} fraction of segments predicted correctly'

    fully_correct_clones = scores_all_tumours.groupby(['clone','tumour_id']).apply(findClonesFullyCorrect).reset_index()
    fully_correct_clones = f'{fully_correct_clones[0].sum() / fully_correct_clones.shape[0]} fraction of clones predicted correctly'

    # labels for illustrator:
    zero = f' Zero predicted as zero {((scores_all_tumours.cn_alpaca == 0) & (scores_all_tumours.cn_true == 0)).sum()}'
    one = f' One predicted as one {((scores_all_tumours.cn_alpaca == 1 )& (scores_all_tumours.cn_true == 1)).sum()}'
    two = f' Two predicted as two {((scores_all_tumours.cn_alpaca == 2) & (scores_all_tumours.cn_true == 2)).sum()}'
    three = f' Three predicted as three {((scores_all_tumours.cn_alpaca == 3) & (scores_all_tumours.cn_true == 3)).sum()}'
    four = f' Four predicted as four {((scores_all_tumours.cn_alpaca == 4) & (scores_all_tumours.cn_true == 4)).sum()}'

    # labels for illustrator - fractions:
    zero_fr = f' Zero predicted as zero {((scores_all_tumours.cn_alpaca == 0) & (scores_all_tumours.cn_true == 0)).sum()/(scores_all_tumours.cn_true == 0).sum()}'
    one_fr = f' One predicted as one {((scores_all_tumours.cn_alpaca == 1 )& (scores_all_tumours.cn_true == 1)).sum()/(scores_all_tumours.cn_true == 1).sum()}'
    two_fr = f' Two predicted as two {((scores_all_tumours.cn_alpaca == 2) & (scores_all_tumours.cn_true == 2)).sum()/(scores_all_tumours.cn_true == 2).sum()}'
    three_fr = f' Three predicted as three {((scores_all_tumours.cn_alpaca == 3) & (scores_all_tumours.cn_true == 3)).sum()/(scores_all_tumours.cn_true == 3).sum()}'
    four_fr = f' Four predicted as four {((scores_all_tumours.cn_alpaca == 4) & (scores_all_tumours.cn_true == 4)).sum()/(scores_all_tumours.cn_true == 4).sum()}'
    with open('sankey_plot_publication_metrics.txt','w') as f:
        f.write('Fully correct segments:'+'\n')
        f.write(fully_correct_segments+'\n')
        f.write('Fully correct clones:'+'\n')
        f.write(fully_correct_clones+'\n')
        f.write('labels for illustrator:'+'\n')
        f.write('zero as zero: '+zero+'\n')
        f.write('one as one: '+one+'\n')
        f.write('two as two: '+two+'\n')
        f.write('three as three: '+three+'\n')
        f.write('four as four: '+four+'\n')
        f.write('labels for illustrator - fractions:'+'\n')
        f.write('zero as zero: '+zero_fr+'\n')
        f.write('one as one: '+one_fr+'\n')
        f.write('two as two: '+two_fr+'\n')
        f.write('three as three: '+three_fr+'\n')
        f.write('four as four: '+four_fr+'\n')

if __name__ == "__main__":
    main()