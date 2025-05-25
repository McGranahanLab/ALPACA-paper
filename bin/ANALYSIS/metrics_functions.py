from scipy.stats import wasserstein_distance as wd
import pandas as pd

def get_total_var_dist(df):
    return 0.5 * (df['proportion_pred'] - df['proportion_true']).abs().sum()


def get_distance_metrics(comparison_df_A,comparison_df_B,tumour_id,segment):
    # calculate total variation distance
    total_var_dist_A = get_total_var_dist(comparison_df_A)
    total_var_dist_B = get_total_var_dist(comparison_df_B)
    total_var_dist_mean = (total_var_dist_A + total_var_dist_B) / 2
    output_df_columns = ['tumour_id','segment','total_var_dist_mean']
    output_df_values = [tumour_id,segment,total_var_dist_mean]
    output_df = pd.DataFrame([output_df_values],columns=output_df_columns)
    return output_df


def get_comparison_dfs(true_copynumbers_seg,true_proportions,model_copynumbers_seg,model_proportions,copy_number_columns=['pred_CN_A','pred_CN_B']):
    all_solutions_seg_complexity_prop = model_copynumbers_seg.merge(model_proportions, on='clone')[copy_number_columns+['proportion']]
    true_copynumbers_seg_prop = true_copynumbers_seg.merge(true_proportions, on='clone')[copy_number_columns+['clone','proportion']]
    output = []
    for copy_number_column in copy_number_columns:
        all_solutions_seg_complexity = all_solutions_seg_complexity_prop[[copy_number_column,'proportion']].groupby(copy_number_column)['proportion'].sum().reset_index()
        true_copynumbers_seg_complexity = true_copynumbers_seg_prop[[copy_number_column,'proportion']].groupby(copy_number_column)['proportion'].sum().reset_index()
        comparison_df = all_solutions_seg_complexity.merge(true_copynumbers_seg_complexity, on=copy_number_column, suffixes=('_pred','_true'),how='outer').fillna(0).rename(columns={copy_number_column:'pred_CN'})
        # make sure proportions sum to 1
        for model in ['pred','true']:
            assert round(comparison_df[f'proportion_{model}'].sum(),3) == 1, f'Proportions do not sum to 1 in allele {copy_number_column}'
        output.append(comparison_df)
    # for backwards compatibility:
    if len(output) == 2:
        return output[0],output[1]
    return output[0]
