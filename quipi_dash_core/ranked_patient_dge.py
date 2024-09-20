import shared as sh

import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go

import numpy as np
import pandas as pd

from statsmodels.stats.multitest import multipletests
from scipy.stats import ranksums

import quipi_dash_core.gene_factor as gf


def fc_coloring(df, fc_thresh, p_val_thresh):
      if df["fc"] < -fc_thresh and df["log10_p_val"] > p_val_thresh:
            return "blue"
      elif df["fc"] > fc_thresh and df["log10_p_val"] > p_val_thresh:
            return "red"
      else:
            return "black"


def factor_score_ranked_dge(gfs_genes, gfs_compartment, quantile, dge_compartment,fc_threshold=1, p_val_thresh=0.000001):
   
   gfs = gf.calculate_gene_factor_score(gfs_genes, gfs_compartment)


   top_patients = gfs[gfs["factor_score"] > gfs["factor_score"].quantile(1-quantile)][["sample_name", "factor_score"]]
   bot_patients = gfs[gfs["factor_score"] < gfs["factor_score"].quantile(quantile)][["sample_name", "factor_score"]]


   top_tpm = sh.quipi_raw[sh.quipi_raw["sample_name"].isin(top_patients["sample_name"])]
   top_comp = top_tpm[top_tpm["compartment"] == dge_compartment]

   bot_tpm = sh.quipi_raw[sh.quipi_raw["sample_name"].isin(bot_patients["sample_name"])]
   bot_comp = bot_tpm[bot_tpm["compartment"] == dge_compartment]

   top_data = top_comp[sh.genes]
   bot_data = bot_comp[sh.genes]

   p_vals = ranksums(top_data,bot_data)[1]
   p_adj = multipletests(p_vals, method = "fdr_bh")[1]

   top_avg_tpm = top_data.mean(axis=0) + .01
   bot_avg_tpm = bot_data.mean(axis=0) + .01

   fc = np.log2(top_avg_tpm / bot_avg_tpm)
   log_p_val = -np.log10(p_adj)

   change_df = pd.DataFrame({"gene":fc.index,
                             "fc":fc,
                             "log10_p_val":log_p_val})
   
   change_df["pt_color"] = change_df.apply(fc_coloring,axis=1,args = (fc_threshold,p_val_thresh))

   fig = px.scatter(change_df,x = "fc", y = "log10_p_val",
      hover_data={"Gene": fc.index},
      labels = {"fc":"Log2(FC)", "log10_p_val": "-log10(P-Value)"},
      color="pt_color",
      color_discrete_map={"red":"red", "black":"black","blue":"blue"})
   
   fig.add_shape(
    type="line",
    x0=fc_threshold, x1=fc_threshold,  # Start and end points of the x-axis
    y0=-.00001, y1=max(log_p_val) * 1.1,  # Cover the entire y-axis range
    line=dict(color="Red", width=2, dash="dash"),  # Line style
    )
   
   fig.add_shape(
    type="line",
    x0=-fc_threshold, x1=-fc_threshold,  # Start and end points of the x-axis
    y0=-.00001, y1=max(log_p_val) * 1.1,  # Cover the entire y-axis range
    line=dict(color="Red", width=2, dash="dash"),  # Line style
    )

   fc_abs_max = abs(max(fc,key=abs))
   fig.update_layout(xaxis=dict(range=[-fc_abs_max-1, fc_abs_max+1]))
   fig.update_layout(showlegend=False)

   fig.add_shape(
    type="line",
    x0=-fc_abs_max-1, x1=fc_abs_max+1,  # Start and end points of the x-axis
    y0=p_val_thresh, y1=p_val_thresh,  # Cover the entire y-axis range
    line=dict(color="Red", width=2, dash="dash"),  # Line style
    )

   fig.update_layout(template='simple_white')

   return fig
        


def feature_ranked_dge(feature_score,compartment,quantile, fc_threshold= 1, p_val_thresh=0.000001 ):
   rank_cat = sh.feature_scores[feature_score]
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
   log_p_val = -np.log10(p_adj)

   change_df = pd.DataFrame({"gene":fc.index,
                             "fc":fc,
                             "log10_p_val":log_p_val})
   
   change_df["pt_color"] = change_df.apply(fc_coloring,axis=1,args = (fc_threshold,p_val_thresh))

   fig = px.scatter(x = fc, 
      y = -np.log10(p_adj),
      hover_data={"Gene": fc.index})
   
   fig = px.scatter(change_df,x = "fc", y = "log10_p_val",
      hover_data={"Gene": fc.index},
      labels = {"fc":"Log2(FC)", "log10_p_val": "-log10(P-Value)"},
      color="pt_color",
      color_discrete_map={"red":"red", "black":"black","blue":"blue"})
   
   fig.add_shape(
    type="line",
    x0=fc_threshold, x1=fc_threshold,  # Start and end points of the x-axis
    y0=-.00001, y1=max(log_p_val) * 1.1,  # Cover the entire y-axis range
    line=dict(color="Red", width=2, dash="dash"),  # Line style
    )
   
   fig.add_shape(
    type="line",
    x0=-fc_threshold, x1=-fc_threshold,  # Start and end points of the x-axis
    y0=-.00001, y1=max(log_p_val) * 1.1,  # Cover the entire y-axis range
    line=dict(color="Red", width=2, dash="dash"),  # Line style
    )

   #fig.update_layout(autosize=False, width=700, height=700)
   fc_abs_max = abs(max(fc,key=abs))
   fig.update_layout(xaxis=dict(range=[-fc_abs_max-1, fc_abs_max+1]))
   fig.update_layout(template='simple_white')
   fig.update_layout(showlegend=False)

   fig.add_shape(
    type="line",
    x0=-fc_abs_max-1, x1=fc_abs_max+1,  # Start and end points of the x-axis
    y0=p_val_thresh, y1=p_val_thresh,  # Cover the entire y-axis range
    line=dict(color="Red", width=2, dash="dash"),  # Line style
    )

   return fig,top_patients,bot_patients
