import pandas as pd
import plotly.express as px
import scanpy as sc
import matplotlib.pyplot as plt
#sc.settings.autoshow = False

def plot_sc_box(gene, x_cat, x_cat_subset, groupby, splitby):
    cols = {gene, x_cat, groupby, splitby} - {"---"}
    input_arr = pd.read_feather("./quipi_humu_data/quipi_humu_adata_clean_full_PROC.feather", columns = cols)
    input_arr = input_arr[input_arr[x_cat].isin(x_cat_subset)]

    splitby = splitby if splitby != "---" else None
    groupby = groupby if groupby != "---" else None

    fig = px.box(input_arr, x = x_cat, y = gene, color = groupby, facet_col=splitby, 
                    facet_col_wrap=3,
                    points = False,
                    facet_col_spacing=0,
                    facet_row_spacing=0,
                    )

    return fig

def plot_sc_dotplot(genes, groupby, splitby, swap):
    adata = sc.read_h5ad("./quipi_humu_data/quipi_humu_adata_clean_full.h5ad", backed="r")

    fig, ax = plt.subplots()

    if groupby == splitby:
        vars = groupby
    else:
        vars = [groupby, splitby] if splitby != "---" else groupby
    
    if swap:
        sc.pl.DotPlot(adata, var_names=genes, groupby=vars, ax=ax).swap_axes().make_figure()
    else:
        sc.pl.DotPlot(adata, var_names=genes, groupby=vars, ax=ax).make_figure()

    return fig

