import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib
import quipi_shared as qsh
import humu_shared as hsh

import box_viol_expression_plot as hbv

def plot_sc_box(gene, x_cat, x_cat_subset, groupby, splitby):
    cols = {gene, x_cat, groupby, splitby} - {"---"}
    input_arr = pd.read_feather("./quipi_humu_data/quipi_humu_adata_clean_full_PROC.feather", columns = cols)
    input_arr = input_arr[input_arr[x_cat].isin(x_cat_subset)]

    splitby = splitby if splitby != "---" else None
    groupby = groupby if groupby != "---" else None

    fig = px.box(input_arr, x = x_cat, y = gene, color = groupby, facet_col=splitby, 
                    facet_col_wrap=2,
                    points = "outliers",
                    facet_col_spacing=0.001,
                    facet_row_spacing=0.03,
                    )

    return fig

def plot_sc_dotplot(genes, groupby, groups, splitby, splits, swap):

    if len(splits) == 0:
        splitby=groupby
        splits = groups

    adata = sc.read_h5ad("./quipi_humu_data/quipi_humu_adata_clean_full.h5ad", backed="r")

    fig, ax = matplotlib.pyplot.subplots()
    if splitby != "---":
        adata = adata[adata.obs[splitby].isin(splits) & (adata.obs[groupby].isin(groups))]
    else:
        adata = adata[adata.obs[groupby].isin(groups)]
    

    if groupby == splitby:
        vars = groupby
    else:
        vars = [groupby, splitby] if splitby != "---" else groupby
    if swap:
        sc.pl.DotPlot(adata, var_names=genes, groupby=vars, ax=ax).swap_axes().make_figure()
    else:
        sc.pl.DotPlot(adata, var_names=genes, groupby=vars, ax=ax).make_figure()

    return fig


def humu_box_comparison_human(human_gene, human_x, human_x_filter, human_groupby, human_transformation, human_comparts):

    fig = hbv.box_viol_exprn(human_transformation, human_x, human_x_filter, human_gene, human_groupby, human_comparts)
    if human_comparts is None:
        fig.update_layout(title= human_transformation + "(" + " ".join(human_gene) + ")")
    else:
        fig.update_layout(title= "Gene-Signature Score(" + " ".join(human_gene) + ") In Compartments: " + " ".join(human_comparts) )

    return fig

def humu_box_comparison_mouse(mouse_gene):
    humu_arr = pd.read_feather("./quipi_humu_data/quipi_humu_adata_clean_full_PROC.feather", columns = hsh.categoricals + [mouse_gene])
    fig = px.box(humu_arr, x = "PanCan_Compartment", y = mouse_gene)
    return fig

    

