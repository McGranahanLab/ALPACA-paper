#!/usr/bin/env python3
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import plotly.express as px


def make_barplot_comparison_plot(alpaca_scores_subclonal,simple_scores_subclonal,alpaca_scores_clonal,simple_scores_clonal):
    font_name = 'DejaVu Sans' # TODO change to Helvetica in local or install Helvetica on HPC
    DPI=300
    plot_size_inches = (3.25, 2.6)
    fig_width, fig_height = plot_size_inches[0]* DPI, plot_size_inches[1]* DPI
    font_size_points = 8
    font_size = (font_size_points/ 72) * DPI
    #subclonal
    alpaca_scores_subclonal = alpaca_scores_subclonal[['tumour_id','clone','clone_accuracy']].copy()
    simple_scores_subclonal = simple_scores_subclonal[['tumour_id','clone','clone_accuracy']].copy()
    alpaca_scores_subclonal.rename(columns={'clone_accuracy':'alpaca_clone_accuracy'}, inplace=True)
    simple_scores_subclonal.rename(columns={'clone_accuracy':'naive_clone_accuracy'}, inplace=True)
    scores_subclonal = pd.merge(alpaca_scores_subclonal, simple_scores_subclonal, on=['tumour_id','clone'])
    scores_subclonal = scores_subclonal[['tumour_id','clone','alpaca_clone_accuracy','naive_clone_accuracy']]
    #clonal
    alpaca_scores_clonal = alpaca_scores_clonal[['tumour_id','clone','clone_accuracy']].copy()
    simple_scores_clonal = simple_scores_clonal[['tumour_id','clone','clone_accuracy']].copy()
    alpaca_scores_clonal.rename(columns={'clone_accuracy':'alpaca_clone_accuracy'}, inplace=True)
    simple_scores_clonal.rename(columns={'clone_accuracy':'naive_clone_accuracy'}, inplace=True)
    scores_clonal = pd.merge(alpaca_scores_clonal, simple_scores_clonal, on=['tumour_id','clone'])
    scores_clonal = scores_clonal[['tumour_id','clone','alpaca_clone_accuracy','naive_clone_accuracy']]
    ########
    alpaca_better_clonal = ((scores_clonal.alpaca_clone_accuracy > scores_clonal.naive_clone_accuracy).sum())/scores_clonal.shape[0]
    naive_better_clonal = ((scores_clonal.alpaca_clone_accuracy < scores_clonal.naive_clone_accuracy).sum())/scores_clonal.shape[0]
    same_score_clonal = ((scores_clonal.alpaca_clone_accuracy == scores_clonal.naive_clone_accuracy).sum())/scores_clonal.shape[0]
    df_clonal = pd.DataFrame({'Accuracy per clone':['Higher with ALPACA','Higher with Simple Model','No difference'],'Fraction of clones':[alpaca_better_clonal,naive_better_clonal,same_score_clonal]})
    df_clonal['Heterogeneity'] = 'Clonal'
    alpaca_better_subclonal = ((scores_subclonal.alpaca_clone_accuracy > scores_subclonal.naive_clone_accuracy).sum())/scores_subclonal.shape[0]
    naive_better_subclonal = ((scores_subclonal.alpaca_clone_accuracy < scores_subclonal.naive_clone_accuracy).sum())/scores_subclonal.shape[0]
    same_score_subclonal = ((scores_subclonal.alpaca_clone_accuracy == scores_subclonal.naive_clone_accuracy).sum())/scores_subclonal.shape[0]
    df_subclonal = pd.DataFrame({'Accuracy per clone':['Higher with ALPACA','Higher with Simple Model','No difference'],'Fraction of clones':[alpaca_better_subclonal,naive_better_subclonal,same_score_subclonal]})
    df_subclonal['Heterogeneity'] = 'Subclonal'
    df = pd.concat([df_clonal,df_subclonal])
    ##########
    df['Fraction of clones'] = df['Fraction of clones'].round(2)
    fig = px.bar(df, x="Heterogeneity", y="Fraction of clones", color="Accuracy per clone", text_auto=True,height=fig_height, width=fig_width,color_discrete_sequence=['#1F78B4', '#FB9A99', '#A8A8A8']) 
    fig.update_layout(font=dict(size=font_size,))
    return fig

def make_barplot_comparison_plot_per_clone(alpaca_scores_per_clone,simple_scores_per_clone,font_size_points=5):
    font_name = 'DejaVu Sans' # TODO change to Helvetica in local or install Helvetica on HPC
    DPI=300
    plot_size_inches = (1, 2)
    fig_width, fig_height = plot_size_inches[0]* DPI, plot_size_inches[1]* DPI
    font_size = (font_size_points/ 72) * DPI
    #subclonal
    alpaca_scores_per_clone = alpaca_scores_per_clone[['tumour_id','clone','clone_accuracy']].copy()
    simple_scores_per_clone = simple_scores_per_clone[['tumour_id','clone','clone_accuracy']].copy()
    alpaca_scores_per_clone.rename(columns={'clone_accuracy':'alpaca_clone_accuracy'}, inplace=True)
    simple_scores_per_clone.rename(columns={'clone_accuracy':'naive_clone_accuracy'}, inplace=True)
    scores_merged = pd.merge(alpaca_scores_per_clone, simple_scores_per_clone, on=['tumour_id','clone'])
    scores_merged = scores_merged[['tumour_id','clone','alpaca_clone_accuracy','naive_clone_accuracy']]
    scores_merged_diff = scores_merged[scores_merged.alpaca_clone_accuracy!=scores_merged.naive_clone_accuracy]
    alpaca_better = ((scores_merged_diff.alpaca_clone_accuracy > scores_merged_diff.naive_clone_accuracy).sum())/scores_merged_diff.shape[0]
    simple_better = ((scores_merged_diff.alpaca_clone_accuracy < scores_merged_diff.naive_clone_accuracy).sum())/scores_merged_diff.shape[0]
    df_plot = pd.DataFrame({'Accuracy per clone':['<br>Higher with<br>ALPACA','Higher with<br>Simple Model'],'Fraction of clones':[alpaca_better,simple_better]})
    df_plot['Fraction of clones'] = df_plot['Fraction of clones'].round(2)
    fig = px.bar(df_plot,x="Accuracy per clone", y="Fraction of clones", color="Accuracy per clone", text_auto=True,height=fig_height, width=fig_width,color_discrete_sequence=['#1F78B4', '#FB9A99']) 
    fig.update_layout(font=dict(size=font_size,))
    fig.update_layout(showlegend=False)
    fig.update_yaxes(tickvals=[0, 0.5, 1])
    fig.update_layout(
    paper_bgcolor='white',
    plot_bgcolor='white')
    return fig



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--alpaca_scores_subclonal', type=str, required=True, help='Alpaca results scord against known copy-numbers')
    parser.add_argument('--alpaca_scores_clonal', type=str, required=True, help='Alpaca results scord against known copy-numbers')
    parser.add_argument('--simple_scores_subclonal', type=str, required=True, help='Simple model results scord against known copy-numbers')
    parser.add_argument('--simple_scores_clonal', type=str, required=True, help='Simple model results scord against known copy-numbers')

    args = parser.parse_args()
    alpaca_scores_subclonal = pd.read_csv(args.alpaca_scores_subclonal)
    alpaca_scores_clonal = pd.read_csv(args.alpaca_scores_clonal)
    simple_scores_subclonal = pd.read_csv(args.simple_scores_subclonal)
    simple_scores_clonal = pd.read_csv(args.simple_scores_clonal)
    fig = make_barplot_comparison_plot(alpaca_scores_subclonal,simple_scores_subclonal,alpaca_scores_clonal,simple_scores_clonal)
    fig.write_image(f'barplot_clonal_subclonal_accuracy.pdf')

if __name__ == "__main__":
    main()