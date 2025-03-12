import pandas as pd
import shared as sh

import plotly.express as px

def plot_sc_violin(gene, x_cat, x_cat_subset, groupby, splitby):
    cols = {gene, x_cat, groupby, splitby} - {"---"}
    input_arr = pd.read_feather("./quipi_humu_data/quipi_humu_adata_clean_full_PROC.feather", columns = cols)
    input_arr = input_arr[input_arr[x_cat].isin(x_cat_subset)]

    splitby = splitby if splitby != "---" else None
    groupby = groupby if groupby != "---" else None

    fig = px.violin(input_arr, x = x_cat, y = gene, color = groupby, facet_col=splitby, facet_col_wrap=3)

    return fig