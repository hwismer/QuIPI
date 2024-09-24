import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go

import numpy as np
import pandas as pd

import shared as sh


def box_viol_exprn(transform, x_cat, gene, groupby, plot_type):

    input_arr = sh.transformations[transform]
    group = sh.categoricals_dict[groupby]

    input_arr["gene_mean"] = input_arr[list(gene)].mean(axis=1)

    if len(gene) != 0:

        if plot_type == "Boxplot":
            fig = px.box(input_arr, x = x_cat, y = "gene_mean", color = group,
                        color_discrete_sequence=px.colors.qualitative.D3,
                        labels=sh.categoricals_dict_reversed)
        else:
            fig = px.violin(input_arr, x = x_cat, y = "gene_mean", color = group,
                            color_discrete_sequence=px.colors.qualitative.D3,
                            labels=sh.categoricals_dict_reversed)
            
        fig.update_layout(title_text= "Mean(" + " ".join(gene) + ") " + transform + "(TPM)", title_x=0.5)
        return fig
