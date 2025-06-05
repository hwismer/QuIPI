import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go

import numpy as np
import pandas as pd
import humu_shared as hsh


def box_humu_flow(score1, score2, x_cat, x_cat_filts, group, sort):
    cols = set([col for col in [score1, score2, x_cat, group] if col in hsh.flow_cats + hsh.flow_scores])
    flow_table = pd.read_feather("./quipi_humu_data/quipi_humu_flow_table.feather", columns=cols)
    flow_table = flow_table[flow_table[x_cat].isin(x_cat_filts)]

    if group in hsh.flow_cats and group != x_cat:
        group = group
    else:
        group = None

    if score2 != "---":
        col_name = "Log2(" + score1 + " / " + score2 + ")"
        flow_table[col_name] = np.log2(flow_table[score1] / flow_table[score2])
    else:
        col_name = score1
    
    if sort is True:
        orders = list(flow_table.groupby(x_cat)[col_name].median().sort_values(ascending=False).index)
    else:
        orders = hsh.humu_flow_orders[x_cat]

    
    fig = px.box(flow_table, x = x_cat , y = col_name, color = group, category_orders={x_cat:orders})



    return fig

