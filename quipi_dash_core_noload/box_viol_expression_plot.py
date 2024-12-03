import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go

import numpy as np
import pandas as pd

import shared as sh


def box_viol_exprn(transform, x_cat, genes, groupby, plot_type):

    if transform == "Raw":
        input_arr = pd.read_csv("./quipi_raw_tpm.csv", usecols=sh.non_genes + genes)
    elif transform == "Log2":
        input_arr = pd.read_csv("./quipi_log2_tpm.csv", usecols=sh.non_genes + genes)

    group = sh.categoricals_dict[groupby]

    input_arr["gene_mean"] = input_arr[list(genes)].mean(axis=1)

    if len(genes) != 0:

        if plot_type == "Boxplot":
            fig = px.box(input_arr, x = x_cat, y = "gene_mean", color = group,
                        color_discrete_sequence=px.colors.qualitative.D3,
                        labels=sh.categoricals_dict_reversed)
        else:
            fig = px.violin(input_arr, x = x_cat, y = "gene_mean", color = group,
                            color_discrete_sequence=px.colors.qualitative.D3,
                            labels=sh.categoricals_dict_reversed)
            
        fig.update_layout(title_text= "Mean(" + " ".join(genes) + ") " + transform + "(TPM)", title_x=0.5)
        return fig
