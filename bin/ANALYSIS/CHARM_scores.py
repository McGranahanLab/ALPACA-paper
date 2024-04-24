#!/usr/bin/env python3
from scipy.stats import pearsonr
import pandas as pd
import plotly.express as px
import argparse
import os

def make_plot(df,xaxis_title,p_value,correlation_coefficient,font_size=22):
    fig = px.scatter(df, x="fraction", y="CharmTSG_score", trendline="ols")
    fig.update_traces(marker=dict(color='black'), selector=dict(mode='markers'))
    fig.update_traces(line=dict(color='black'), selector=dict(mode='lines'))
    fig.update_layout(
        plot_bgcolor='white',  
        paper_bgcolor='white',
        xaxis_title=xaxis_title,
        yaxis_title='CharmTSG Score',  
        font=dict(
            family='Helvetica',  
            size=font_size,
            color='black',
        )
    )
    fig.add_annotation(
        x=0.05, y=0.95,
        xref="paper",
        yref="paper",  
        text=f"P-value: {round(p_value,3)}",
        showarrow=False,
        font=dict(
            family="Helvetica",
            size=font_size-2,
            color='black',
        )
    )
    fig.add_annotation(
        x=0.05, y=0.9,
        xref="paper",
        yref="paper",  
        text=f"R: {round(correlation_coefficient,3)}",
        showarrow=False,
        font=dict(
            family="Helvetica",
            size=font_size-2,
            color='black',
        )
    )
    for i, row in df.iterrows():
        fig.add_annotation(
            x=row['fraction'],
            y=row['CharmTSG_score'],
            text=row['chr_arm'],
            showarrow=False,
            font=dict(
                family="Helvetica",
                size=font_size-4,
            ),
            yshift=10
        )
    fig.update_layout(
    xaxis_showgrid=True, 
    yaxis_showgrid=True,  
    xaxis_gridcolor='Black',  
    yaxis_gridcolor='Black',  
    )
    return fig


def process_charm_table(df_path):
    # load Davoli table:
    # Supplementary table 6, tab Table S6B
    # https://www.sciencedirect.com/science/article/pii/S0092867413012877
    df = pd.read_excel(df_path,sheet_name='Table S6B', header=2)
    df = df[['Arm','Density_TSG_in_the_arm','Average_Deletion_Frequency','CharmTSG_score']]
    return df


def get_correlation_and_P(df):
    # correlate with CHARM score:
    correlation_coefficient, p_value = pearsonr(df['CharmTSG_score'], df['fraction'])
    print(f"Correlation Coefficient: {correlation_coefficient}")
    print(f"P-value: {p_value}")
    return correlation_coefficient, p_value


def find_regional_arm_level_loh(input_df):
    input_df['A_lost'] = sum((input_df['cpnA']<0.5)*input_df['seg_len'])/(input_df['seg_len'].sum())
    input_df['B_lost'] = sum((input_df['cpnB']<0.5)*input_df['seg_len'])/(input_df['seg_len'].sum())
    return input_df


def make_plotting_df(regional_arm_level_events_loh,charm_table,colname):
    arms_only_alpaca_unique = regional_arm_level_events_loh[['tumour_id','chr_arm',colname]].drop_duplicates()
    arms_only_alpaca_unique.rename(columns={colname:'has_arm_loh'},inplace=True)
    total_arm_event_count = pd.DataFrame(arms_only_alpaca_unique.groupby('chr_arm')['has_arm_loh'].sum())
    total_arm_event_count['fraction'] = total_arm_event_count['has_arm_loh'] / len(arms_only_alpaca_unique.tumour_id.unique())
    total_arm_event_count.reset_index(inplace=True)
    df = pd.merge(total_arm_event_count[['chr_arm','fraction']],charm_table,left_on='chr_arm',right_on='Arm')
    df['fraction'] = df['fraction'].astype(float)
    df['CharmTSG_score'] = df['CharmTSG_score'].astype(float)
    return df


def make_input_table(input_data_directory,combined_results,arm_level_events_alpaca):
    fractional_input = pd.concat([pd.read_csv(f'{input_data_directory}/{tumour_id}/ALPACA_input_table.csv') for tumour_id in os.listdir(input_data_directory) if 'CRUK' in tumour_id])
    fractional_input = fractional_input[fractional_input.segment.isin(combined_results.segment.unique())]
    arm_level_events_alpaca = arm_level_events_alpaca[['tumour_id','segment','chr_arm','any_arm_loh_acquired_edge']].drop_duplicates()
    combined_results = combined_results[['tumour_id','segment','pred_CN_A','pred_CN_B']].drop_duplicates()
    alpaca_with_regional = pd.merge(combined_results,fractional_input,left_on=['tumour_id','segment'],right_on=['tumour_id','segment'])
    regional_arm_level_events = pd.merge(alpaca_with_regional,arm_level_events_alpaca,left_on=['tumour_id','segment'],right_on=['tumour_id','segment'])
    regional_arm_level_events['seg_len'] = regional_arm_level_events.segment.str.split('_').str[2].astype(int)-regional_arm_level_events.segment.str.split('_').str[1].astype(int)
    regional_arm_level_events_loh = regional_arm_level_events.groupby(['tumour_id','sample','chr_arm']).apply(find_regional_arm_level_loh)
    regional_arm_level_loss_thr = 0.98
    regional_arm_level_events_loh['regional_arm_level_loh'] = (regional_arm_level_events_loh.A_lost>regional_arm_level_loss_thr)|(regional_arm_level_events_loh.B_lost>regional_arm_level_loss_thr)
    regional_arm_level_events_loh.rename(columns={'any_arm_loh_acquired_edge':'alpaca_arm_level_loh'},inplace=True)
    regional_arm_level_events_loh['alpaca_unique'] = (regional_arm_level_events_loh.alpaca_arm_level_loh)&(~regional_arm_level_events_loh.regional_arm_level_loh)
    regional_arm_level_events_loh['region_unique'] = (~regional_arm_level_events_loh.alpaca_arm_level_loh)&(regional_arm_level_events_loh.regional_arm_level_loh)
    return regional_arm_level_events_loh

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm_level_events_alpaca",type=str,required=True,help='')
    parser.add_argument("--combined_results", type=str, required=True, help='')
    parser.add_argument("--input_data_directory",type=str,required=True,help='')
    parser.add_argument("--charm_table", type=str, required=True, help='')
    args = parser.parse_args()

    arm_level_events_alpaca = pd.read_csv(args.arm_level_events_alpaca)
    combined_results = pd.read_csv(args.combined_results)
    input_data_directory = args.input_data_directory
    
    charm_table = process_charm_table(args.charm_table)
    regional_arm_level_events_loh = make_input_table(input_data_directory,combined_results,arm_level_events_alpaca)
    
    alpaca_unique = make_plotting_df(regional_arm_level_events_loh.copy(),charm_table.copy(),'alpaca_unique') # only ALPACA unique calls, i.e. not visible regionally
    region_unique = make_plotting_df(regional_arm_level_events_loh.copy(),charm_table.copy(),'regional_arm_level_loh') # all regional calls, regardless of the detection by ALPACA
    alpaca_all = make_plotting_df(regional_arm_level_events_loh.copy(),charm_table.copy(),'alpaca_arm_level_loh') # all ALAPCA calls including those that are also visible regionally

    alpaca_unique_correlation_coefficient, alpaca_unique_p_value = get_correlation_and_P(alpaca_unique)
    region_unique_correlation_coefficient, region_unique_p_value = get_correlation_and_P(region_unique)
    alpaca_all_correlation_coefficient, alpaca_all_p_value = get_correlation_and_P(alpaca_all)

    alpaca_unique_fig = make_plot(alpaca_unique,'Fraction of arm-level LOH detected with ALPACA',alpaca_unique_p_value,alpaca_unique_correlation_coefficient,30)
    alpaca_unique_fig.write_image('correlation_with_charm_tsg_score_unique_ALPACA_LOH_arm_level.png',width=2000, height=800, scale=2)

    region_unique_fig = make_plot(region_unique,'Fraction of arm-level LOH detected regionally',region_unique_p_value,region_unique_correlation_coefficient,30)
    region_unique_fig.write_image('correlation_with_charm_tsg_score_region_LOH_arm_level.png',width=2000, height=800, scale=2)

    alpaca_all_fig = make_plot(alpaca_all,'Fraction of arm-level LOH detected with ALPACA',alpaca_all_p_value,alpaca_all_correlation_coefficient,30)
    alpaca_all_fig.write_image('correlation_with_charm_tsg_score_all_ALPACA_LOH_arm_level.png',width=2000, height=800, scale=2)
    
    
#TODO remove dev:
'''
arm_level_events_alpaca = pd.read_csv('/camp/home/pawlikp/CN-CCF/publication/output/primary/primary_default/cohort_outputs/edge_cn_change_seg_chrarm.csv')
combined_results = pd.read_csv('/camp/home/pawlikp/CN-CCF/publication/output/primary/primary_default/cohort_outputs/combined.csv')
input_data_directory = '/camp/home/pawlikp/CN-CCF/publication/input/primary'
charm_table = process_charm_table('/camp/home/pawlikp/CN-CCF/publication/_assets/1-s2.0-S0092867413012877-mmc6.xlsx')
os.chdir('/camp/project/proj-tracerx-lung/tctProjects/CN-CCF/publication/output/primary/primary_default/cohort_outputs/CHARM_scores')
'''






