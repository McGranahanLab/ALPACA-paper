import sys
from functions import pval_to_annotations
from scipy.stats import ttest_rel, wilcoxon, rankdata
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

publication_dir = f"../.."


def get_hd_hatchet_alpaca_colours(palette):
    import json

    with open(palette, "r") as f:
        palette = json.load(f)
    catgorical_palette = palette["categorical"]
    colours = {}
    colours["hd"] = catgorical_palette["c12"]
    colours["hatchet"] = catgorical_palette["c8"]
    colours["alpaca"] = catgorical_palette["c2"]
    return colours


def effect_size_rank_biserial_and_p_val(x, y):
    # Perform Wilcoxon test
    w, p = wilcoxon(x, y)
    # Rank the data
    ranks = rankdata(np.concatenate((x, y)))
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
    ranks = rankdata(np.abs(differences))
    positive_rank_sum = np.sum(ranks[differences > 0])
    negative_rank_sum = np.sum(ranks[differences < 0])
    # Use the smaller of the two rank sums
    R_smaller = min(positive_rank_sum, negative_rank_sum)
    r = 1 - (R_smaller / (n1 * (n1 + 1) / 2))
    return r, p, rank_biserial


def make_error_comparison_plot(
    cohort_results,
    y_axis_title,
    palette_path="../../_assets/publication_palette.json",
    w=600,
    h=800,
    font_size=45,
    orientation="v",
    fonts="Arial",
    test="wilcoxon",
):
    colours = get_hd_hatchet_alpaca_colours(palette_path)
    # wilcoxon test:
    effect_size_hatchet, p_value_wilcoxon, rank_biserial_hatchet = (
        effect_size_rank_biserial_and_p_val(
            cohort_results.alpaca, cohort_results.hatchet
        )
    )
    results = pd.DataFrame(
        {
            "method": ["HATCHet", "ALPACA"],
            "effect_size": [effect_size_hatchet, ""],
            "p_value_wilcoxon": [p_value_wilcoxon, ""],
            "rank_biserial": [rank_biserial_hatchet, ""],
            "mean": [np.mean(cohort_results.hatchet), np.mean(cohort_results.alpaca)],
            "median": [
                np.median(cohort_results.hatchet),
                np.median(cohort_results.alpaca),
            ],
        }
    )

    selected_p_value = p_value_wilcoxon
    fig = go.Figure()

    jitter_value = 0
    point_size = font_size / 10
    line_width = 2
    alpaca_colour = colours["alpaca"]
    hatchet_colour = colours["hatchet"]
    if orientation == "v":
        pval_y = 1.1
        fig.update_xaxes(showgrid=True)
        fig.update_yaxes(
            showgrid=True,
            gridwidth=1,
            gridcolor="black",
            zeroline=True,
            zerolinecolor="black",
            zerolinewidth=1,
            title=y_axis_title,
            title_font=dict(family=fonts, size=font_size * 0.8),
        )

        hatchet_trace = go.Box(
            y=cohort_results.hatchet,
            orientation=orientation,
            boxpoints="outliers",
            jitter=jitter_value,
            pointpos=0,
            line=dict(color=hatchet_colour, width=line_width),
            marker=dict(size=point_size),
            name="HATCHet",
        )
        alpaca_trace = go.Box(
            y=cohort_results.alpaca,
            orientation=orientation,
            boxpoints="outliers",
            jitter=jitter_value,
            pointpos=0,
            line=dict(color=alpaca_colour, width=line_width),
            marker=dict(size=point_size),
            name="ALPACA",
        )
    else:
        pval_y = 1.2
        fig.update_yaxes(showgrid=False)
        fig.update_xaxes(
            showgrid=False,
            gridwidth=1,
            gridcolor="black",
            zeroline=True,
            zerolinecolor="black",
            zerolinewidth=1,
            title=y_axis_title,
            title_font=dict(family=fonts, size=font_size * 0.8),
        )

        hatchet_trace = go.Box(
            x=cohort_results.hatchet,
            orientation=orientation,
            boxpoints="outliers",
            jitter=jitter_value,
            pointpos=0,
            line=dict(color=hatchet_colour, width=line_width),
            marker=dict(size=point_size),
            name="HATCHet",
        )
        alpaca_trace = go.Box(
            x=cohort_results.alpaca,
            orientation=orientation,
            boxpoints="outliers",
            jitter=jitter_value,
            pointpos=0,
            line=dict(color=alpaca_colour, width=line_width),
            marker=dict(size=point_size),
            name="ALPACA",
        )
    fig.update_layout(
        title="",
        font=dict(family=fonts, size=font_size),
        width=w,
        height=h,
        showlegend=False,
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        annotations=[
            {
                "text": f"{pval_to_annotations(selected_p_value)}",
                "x": 0.5,
                "y": pval_y,
                "xref": "paper",
                "yref": "paper",
                "align": "left",
                "showarrow": False,
                "font": {"size": 24},
            },
        ],
    )

    fig.add_trace(hatchet_trace)
    fig.add_trace(alpaca_trace)
    fig.update_layout(
        font=dict(family=fonts, size=font_size),
        width=w,
        height=h,
        paper_bgcolor="white",
        plot_bgcolor="white",
    )
    return fig, results


def get_colourmap(max_val, testing=False, cn_type="total"):
    # modified version of get_colourmap found in publication/bin/ANALYSIS/CLUSTERING_make_clustermap_plot.py
    min_val = 0
    positive_colors = np.array([1, 0, 0, 1])
    negative_colors = np.array([0, 0, 1, 1])
    white = np.array([1, 1, 1, 1])
    if cn_type == "total":
        colors = [
            (0.0, negative_colors),  # Blue at 0
            (1.0 / max_val, negative_colors),  # Blue at 1
            (2.0 / max_val, white),  # White at 2
            (1.0, positive_colors),  # Red at max_val
        ]
    else:
        colors = [
            (0.0, negative_colors),  # Blue at 0
            (1.0 / max_val, white),  # White at 2
            (1.0, positive_colors),  # Red at max_val
        ]
    custom_colormap = LinearSegmentedColormap.from_list("custom_map", colors)
    # testing:
    if testing:
        data = np.random.uniform(low=min_val, high=max_val, size=(10, 10))
        data = np.array(range(min_val, max_val + 1)).reshape(1, -1)
        plt.imshow(data, cmap=custom_colormap)
        plt.colorbar()
        plt.show()
    return custom_colormap
