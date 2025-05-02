import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib
import quipi_shared as qsh
import humu_shared as hsh

def plot_sc_box(gene, x_cat, x_cat_subset, groupby, splitby):
    cols = {gene, x_cat, groupby, splitby} - {"---"}
    input_arr = pd.read_feather("./quipi_humu_data/quipi_humu_adata_clean_full_PROC.feather", columns = cols)
    input_arr = input_arr[input_arr[x_cat].isin(x_cat_subset)]

    splitby = splitby if splitby != "---" else None
    groupby = groupby if groupby != "---" else None

    fig = px.violin(input_arr, x = x_cat, y = gene, color = groupby, facet_col=splitby, 
                    facet_col_wrap=2,
                    #violinmode="overlay",
                    #points = "all",
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

def humu_box_comparison(human_gene, mouse_gene):
    fig = make_subplots(rows=1, cols=2,
                        subplot_titles=('Human Log2(TPM) ' + human_gene , 'Murine Log-Normalized Counts ' + mouse_gene))
    print(human_gene, mouse_gene)

    human_arr = pd.read_feather("./quipi_data/quipi_log2_tpm.feather", columns=qsh.non_genes + [human_gene])
    humu_arr = pd.read_feather("./quipi_humu_data/quipi_humu_adata_clean_full_PROC.feather", columns = hsh.categoricals + [mouse_gene])

    ax1 = px.box(human_arr, x = "compartment", y = human_gene)
    ax2 = px.box(humu_arr, x = "PanCan_Compartment", y = mouse_gene)
    fig.add_trace(px.box(human_arr, x = "compartment", y = human_gene,).data[0], row=1,col=1)
    fig.add_trace(px.box(humu_arr, x = "PanCan_Compartment", y = mouse_gene).data[0], row=1, col=2)

    return fig

    

