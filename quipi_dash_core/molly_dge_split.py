import shared as sh

from shiny import App, render, ui, reactive
import plotly.express as px
from shinywidgets import output_widget, render_widget 
from plotly.subplots import make_subplots
import plotly.graph_objects as go

import numpy as np
import pandas as pd

from scipy.stats import zscore
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from scipy.stats import ranksums



def molly_dge(feature_score,compartment,quantile, fc_threshold= 1, p_val_thresh=0.000001 ):
        rank_cat = sh.flow_scores[feature_score]
        comp = compartment
        quantile = quantile

        top_patients = sh.quipi_flow[sh.quipi_flow[rank_cat] > sh.quipi_flow[rank_cat].quantile(1 - quantile)][["sample_name", rank_cat]]
        bot_patients = sh.quipi_flow[sh.quipi_flow[rank_cat] < sh.quipi_flow[rank_cat].quantile(quantile)][["sample_name", rank_cat]]


        top_tpm = sh.quipi_raw[sh.quipi_raw["sample_name"].isin(top_patients["sample_name"])]
        top_comp = top_tpm[top_tpm["compartment"] == comp]

        bot_tpm = sh.quipi_raw[sh.quipi_raw["sample_name"].isin(bot_patients["sample_name"])]
        bot_comp = bot_tpm[bot_tpm["compartment"] == comp]

        top_data = top_comp[sh.genes]
        bot_data = bot_comp[sh.genes]

        p_vals = ranksums(top_data,bot_data)[1]
        p_adj = multipletests(p_vals, method = "fdr_bh")[1]

        top_avg_tpm = top_data.mean(axis=0) + .01
        bot_avg_tpm = bot_data.mean(axis=0) + .01

        fc = np.log2(top_avg_tpm / bot_avg_tpm)

        fig = px.scatter(x = fc, 
           y = -np.log10(p_adj),
           hover_data={"Gene": fc.index})

        #fig.update_layout(autosize=False, width=700, height=700)
        fc_abs_max = abs(max(fc,key=abs))
        fig.update_layout(xaxis=dict(range=[-fc_abs_max-1, fc_abs_max+1]))
        fig.update_layout(template='simple_white')

        return fig,top_patients,bot_patients
