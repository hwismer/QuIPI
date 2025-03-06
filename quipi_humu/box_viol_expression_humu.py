import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go

import numpy as np
import pandas as pd
import shared as sh


def box_humu(score1, score2, x_cat, x_cat_filts, group):
    cols = set([col for col in [score1, score2, x_cat, group] if col in sh.flow_cats + sh.flow_scores])
    flow_table = pd.read_feather("./quipi_humu/quipi_humu_data/quipi_humu_flow_table.feather", columns=cols)
    flow_table = flow_table[flow_table[x_cat].isin(x_cat_filts)]

    if group in sh.flow_cats and group != x_cat:
        group = group
    else:
        group = None

    if score2 != "---":
        flow_table[score1 + " / " + score2] = flow_table[score1] / flow_table[score2]
        fig = px.box(flow_table, x = x_cat , y = score1 + " / " + score2, color = group)
    else:
        print(flow_table, x_cat, score1, group)
        fig = px.box(flow_table, x = x_cat, y = score1, color = group )

    return fig