import json
import os
from datetime import datetime
import pandas as pd
import numpy as np
import math
import kneed
from scipy.stats import norm
import typing
from ALPACA_model_class import Model

def ensure_elbow_strictly_decreasing(df):
    for col in ['D_score']:
        df[col] = df[col].round(3)
        current_minimum = 1000
        new_values = []
        for v in df[col]:
            if v < current_minimum:
                new_values.append(v)
                current_minimum = v
            else:
                new_v = current_minimum - 0.001
                current_minimum = new_v
                new_values.append(new_v)
        df[col] = new_values
    return df


def find_s_values(elbow_search_df, max_iterations, x='allowed_complexity', y='D_score'):
    elbow_table = (elbow_search_df[[x, y]]).dropna()
    m = 1
    s = kneed.KneeLocator(elbow_table[x], elbow_table[y], S=m, curve="convex", direction="decreasing", interp_method="interp1d", online=True).knee
    s_code = 'default'
    if not s:
        # try with higher sensitivity:
        s_candidates = list(
            set([s for s in [kneed.KneeLocator(elbow_table[x], elbow_table[y], S=m, curve="convex", direction="decreasing", interp_method="interp1d", online=True).knee for m in range(0, 200)] if
                 s is not None]))
        if len(s_candidates) > 0:
            s = s_candidates[0]
            s_code = 'high_sensitivity'
    return s, s_code


def find_elbow(list_of_solutions, max_iterations):
    solutions_combined = pd.concat(list_of_solutions)
    elbow_search_df = solutions_combined[['complexity', 'D_score', 'CI_score', 'allowed_complexity']].drop_duplicates(subset='allowed_complexity', keep='first').sort_values(
        'allowed_complexity').reset_index(drop=True)
    elbow_search_df_strictly_decreasing = ensure_elbow_strictly_decreasing(elbow_search_df.copy())
    s_raw, raw_code = find_s_values(elbow_search_df, max_iterations, 'allowed_complexity', 'D_score')
    s_strictly_decreasing, dec_code = find_s_values(elbow_search_df_strictly_decreasing, max_iterations, 'allowed_complexity', 'D_score')
    s_min = min(s_raw, s_strictly_decreasing)
    s_values = {'s_min': s_min, 's_raw': s_raw, 's_strictly_decreasing': s_strictly_decreasing, 'raw_code': raw_code, 'dec_code': dec_code}
    return s_values


def missing_clones_inherit_from_children(optimal_solution, tree, cp_table):
    """
    Only applies when a missing clone has exactly one child., Otherwise parsimony should be enough to find the correct solution.
    """
    missing_clones = list(cp_table.loc[cp_table.sum(1) == 0].index)
    # exclude leaves as they cannot have any children:
    leaves = [x[-1] for x in tree]
    missing_clones = [x for x in missing_clones if x not in leaves]
    if len(missing_clones) > 0:
        for clone in missing_clones:
            tree_branch = [branch for branch in tree if clone in branch][0]
            missing_in_branch = [c for c in tree_branch if c in missing_clones]
            missing_in_branch.reverse()
            for missing_clone in missing_in_branch:
                non_missing_children = []
                for branch in tree:
                    if missing_clone in branch:
                        child_index = branch.index(missing_clone) + 1
                        child = branch[child_index]
                        if child not in non_missing_children:
                            non_missing_children.append(child)
                if len(non_missing_children) == 1:
                    non_missing_child = non_missing_children[0]
                    optimal_solution.loc[optimal_solution.clone == missing_clone, 'pred_CN_A'] = optimal_solution.loc[optimal_solution.clone == non_missing_child, 'pred_CN_A'].values[0]
                    optimal_solution.loc[optimal_solution.clone == missing_clone, 'pred_CN_B'] = optimal_solution.loc[optimal_solution.clone == non_missing_child, 'pred_CN_B'].values[0]
    return optimal_solution


def remove_small_clones(cp_table, tree):
    tree_levels = max([len(branch) for branch in tree])
    MRCA = tree[0][0]
    for level in reversed(range(0, tree_levels)):
        for branch in tree:
            try:
                clone = branch[level]
                for region in cp_table.columns:
                    clone_cp = cp_table.loc[clone, region]
                    if (clone_cp > 0) & (clone_cp < 0.1) & (clone != MRCA):
                        parent = branch[branch.index(clone) - 1]
                        cp_table.loc[parent, region] = cp_table.loc[parent, region] + clone_cp
                        cp_table.loc[clone, region] = 0
            except IndexError:
                pass
    return cp_table


def read_tree_json(json_path: str):
    with open(json_path, 'r') as f:
        tree = json.load(f)
    return tree


def split_input_file_name(input_file_name: str):
    stripped_name = input_file_name.split('ALPACA_input_table_')[1]
    assert stripped_name.count('_') == 3, 'Input name has to many underscores'
    stripped_name = stripped_name.replace('.csv', '')
    t_id = stripped_name.split('_')[0]
    s_name = '_'.join(stripped_name.split('_')[1:])
    return t_id, s_name


def calculate_CI(df_seg_reg, CI):
    for allele in ['A', 'B']:
        data = np.array(df_seg_reg[f'ph_cpn{allele}_vec'])
        # drop nans:
        data = data[~np.isnan(data)]
        m = np.mean(data)
        std = np.std(data, ddof=1)
        lower_CI = norm.ppf(CI / 2, loc=m, scale=std)
        upper_CI = norm.ppf(1 - CI / 2, loc=m, scale=std)
        df_seg_reg[f'lower_CI_{allele}'] = lower_CI
        df_seg_reg[f'upper_CI_{allele}'] = upper_CI
    return df_seg_reg


def get_ci_table(input_table, tumour_dir, segment, CI=0.5):
    # try to read SNP table, if not present, create artificial CI table:
    try:
        asas_table = pd.read_csv(f'{tumour_dir}/asas_table.csv')
        asas_table = asas_table[asas_table.segment == segment]
        ci_table = asas_table.groupby(['sample', 'segment']).apply(lambda df_seg_reg: calculate_CI(df_seg_reg, CI))
        ci_table = ci_table[['sample', 'segment', f'lower_CI_A', f'upper_CI_A', f'lower_CI_B', f'upper_CI_B']].drop_duplicates()
    except FileNotFoundError:
        # TODO devise more sophisticated way to calculate CI
        print('No SNP table found, creating artificial CI table')
        ci_table = input_table[['sample', 'segment']].drop_duplicates().copy()
        for x in ['lower_CI_A', 'upper_CI_A', 'lower_CI_B', 'upper_CI_B']:
            ci_table[x] = 0
        for s in ci_table['sample'].unique():
            A = input_table[input_table['sample'] == s].ph_cpnA_vec.median()
            B = input_table[input_table['sample'] == s].ph_cpnB_vec.median()
            ci_table.loc[ci_table['sample'] == s, f'lower_CI_A'] = A - 0.5
            ci_table.loc[ci_table['sample'] == s, f'lower_CI_B'] = B - 0.5
            ci_table.loc[ci_table['sample'] == s, f'upper_CI_A'] = A + 0.5
            ci_table.loc[ci_table['sample'] == s, f'upper_CI_B'] = B + 0.5
    for allele in ['A', 'B']:
        ci_table[f'lower_CI_{allele}'] = ci_table[f'lower_CI_{allele}'].apply(lambda x: max(x, 0))
        ci_table[f'upper_CI_{allele}'] = ci_table[f'upper_CI_{allele}'].apply(lambda x: max(x, 0.01))
    ci_table['ci_value'] = CI
    ci_table = ci_table.reset_index(drop=True).sort_values('sample')
    return ci_table


def validate_inputs(it: pd.DataFrame, cpt: pd.DataFrame, cit: pd.DataFrame, t: typing.List[typing.List]):
    # check if all clones are present:
    cpt_clones = set(cpt.index.unique())
    tree_clones = set([c for branch in t for c in branch])
    if cpt_clones != tree_clones:
        raise ValueError('Clones in cp_table and tree_paths.json do not match')
    # check if all samples are present:
    it_samples = set(it['sample'].unique())
    cpt_samples = set(cpt.columns)
    cit_samples = set(cit['sample'].unique())
    if it_samples != cpt_samples != cit_samples:
        raise ValueError('Samples in input table, cp_table and ci_table do not match')
    # check if segment is present in the ci_table:
    it_segments = set(it['segment'].unique())
    cit_segments = set(cit['segment'].unique())
    if not it_segments.issubset(cit_segments):
        raise ValueError('Segments in input table and ci_table do not match')


def calibrate_clone_proportions(cp: pd.DataFrame):
    for r in cp.select_dtypes(include=float).columns:
        cp[r] = cp[r] / cp[r].sum()
    return cp


class SegmentSolution:
    def __init__(self, input_file_name, config=None):
        if config is None:
            config = {'preprocessing_config': {}, 'model_config': {}}
        self.ci = 0.5
        self.rsc = False
        self.input_data_directory = None
        self.ccp = True
        self.s_type = 's_strictly_decreasing'
        self.optimal_solution = None
        self.optimal_solution_index = None
        self.elbow = None
        self.maximum_complexity = None
        self.solutions_combined = None
        self.config = config
        self.metrics = None
        self.no_change_in_complexity = False
        self.no_change_in_D_score = False
        self.no_improvement_in_D_score = False
        self.diploid_solution_found = False
        self.compare_with_true_solution = False
        # load config
        for key, value in self.config['preprocessing_config'].items():
            setattr(self, key, value)
        self.tumour_id, self.segment = split_input_file_name(input_file_name)
        # define tumour input directory depending on the run environment:
        if 'NXF_WORK' in os.environ:  # running on Nextflow TODO move out of the class
            self.tumour_dir = f'{self.input_data_directory}/{self.tumour_id}'
        else:
            self.tumour_dir = f'..'
        # load fractional copy numbers:
        self.input_table = pd.read_csv(input_file_name).sort_values('sample')
        # load tree:
        self.tree = read_tree_json(f'{self.tumour_dir}/tree_paths.json')
        # load clone proportions:
        self.cp_table = pd.read_csv(f'{self.tumour_dir}/cp_table.csv', index_col='clone')
        self.cp_table = calibrate_clone_proportions(self.cp_table) if self.ccp else self.cp_table
        self.cp_table = remove_small_clones(self.cp_table, self.tree) if self.rsc else self.cp_table
        # get confidence intervals for copy number values:
        self.ci_table = get_ci_table(self.input_table, self.tumour_dir, self.segment, CI=self.ci)
        # check if inputs are in the expected format and contain all required columns:
        validate_inputs(it=self.input_table, cpt=self.cp_table, cit=self.ci_table, t=self.tree)
        #
        print(datetime.now())
        print(f'Running: {input_file_name}')
        print(f'Tumour id: {self.tumour_id}')
        print(f'Segment name: {self.segment}')
        print('input table:')
        print(self.input_table)
    
    def get_model_metrics(self, model_iteration):
        self.metrics['D_scores'].append(model_iteration.solution.D_score.iloc[0])
        self.metrics['solutions'].append(model_iteration.solution)
        self.metrics['run_time'].append(model_iteration.model.Runtime)
        self.metrics['models'].append(model_iteration.model)
        self.metrics['complexity'].append(model_iteration.solution.complexity.iloc[0])
    
    def run_model(self, allowed_complexity):
        allowed_complexity = {'allowed_tree_complexity': allowed_complexity}
        model_iteration = Model(
            segment=self.segment,
            ci_table=self.ci_table,
            fractional_copy_number_table=self.input_table,
            tree=self.tree,
            clone_proportions=self.cp_table,
            **{**self.config['model_config'], **allowed_complexity})
        model_iteration.model.optimize()
        model_iteration.get_output()
        self.get_model_metrics(model_iteration)
    
    def stop_conditions_check(self, oft):
        optimization_time = self.metrics['run_time'][-1]
        slow_iteration = optimization_time >= self.metrics['models'][-1].params.TimeLimit
        no_improvement_in_D_score = sum(abs(np.diff(self.metrics['D_scores'])) < oft) > 3
        at_least_3_complexities = len(set(self.metrics['complexity'])) >= 3
        return no_improvement_in_D_score & at_least_3_complexities & slow_iteration
    
    def run_iterations(self):
        # use heuristics to determine maximum complexity:
        self.maximum_complexity = max(20, len(self.input_table['sample'].unique()) * math.ceil(self.input_table[['cpnA', 'cpnB']].max().max()))
        objective_function_threshold = 0.1  # iterations will stop if D score does not improve by more than this value in 3 consecutive iterations
        self.metrics = {name: [] for name in ["D_scores", "solutions", "run_time", "models", 'complexity']}
        # run diploid model:
        self.run_model(allowed_complexity=0)
        # don't iterate if solution is likely to be diploid:
        if self.metrics['D_scores'][0] > objective_function_threshold:
            complexity_range = range(1, self.maximum_complexity)
            for c in complexity_range:
                print(f'**Iterating with complexity: {c}')
                self.run_model(allowed_complexity=c)
                stop_conditions = self.stop_conditions_check(objective_function_threshold)
                if stop_conditions:
                    # check if elbow can be found
                    elbow_findable = find_elbow(self.metrics['solutions'], self.maximum_complexity)['s_min'] < 1000  # TODO convert to method
                    if elbow_findable:
                        print(f'** Stopping iterations at complexity {c} due to lack of improvement in D score')
                        break
        self.solutions_combined = pd.concat(self.metrics['solutions'])
    
    def find_optimal_solution(self):
        # check if diploid solution was found:
        diploid_solution_found = len(set(self.metrics['D_scores'])) == 1
        if diploid_solution_found:
            self.optimal_solution_index = 0
        else:
            # find elbow:
            self.elbow = find_elbow(self.metrics['solutions'], self.maximum_complexity)
            self.optimal_solution_index = self.elbow[self.s_type]
    
    def get_solution(self, s=None):
        if s is None:
            s = self.optimal_solution_index
        self.optimal_solution = self.solutions_combined.query(f"allowed_complexity == {s}").copy()
        self.optimal_solution = missing_clones_inherit_from_children(self.optimal_solution, self.tree, self.cp_table)
        self.optimal_solution['tumour_id'] = self.tumour_id
        self.optimal_solution['segment'] = self.segment
        self.optimal_solution.drop(columns=['allowed_complexity'], inplace=True)
