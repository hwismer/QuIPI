from shiny import App, render, ui, reactive
import plotly.express as px
from shinywidgets import output_widget, render_widget 
from plotly.subplots import make_subplots
import plotly.graph_objects as go

import shared as sh

import numpy as np
import pandas as pd
pd.options.plotting.backend = 'plotly'
import matplotlib.pyplot as plt


from scipy.stats import zscore
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from scipy.stats import ranksums

import ranked_patient_dge as rpd
import gene_factor as gf
import box_viol_expression_plot as bv

def plot_pancan_exprn_subplots(transform, genes, compartment):

    #input_arr = sh.transformations[transform]
    if transform == "Raw":
        input_arr = pd.read_csv("./data/quipi_raw_tpm.csv", usecols=sh.non_genes + genes)
    elif transform == "Log2":
        input_arr = pd.read_csv("./data/quipi_log2_tpm.csv", usecols=sh.non_genes + genes)
    input_arr = input_arr[input_arr["compartment"] == compartment]



    if len(genes) != 0:

        n_col = min(4,len(genes))
        n_rows = (len(genes) + n_col - 1) // n_col

        fig = make_subplots(rows =  n_rows , cols = n_col, 
                            subplot_titles=genes,
                            vertical_spacing=.05,horizontal_spacing=.02,
                            shared_xaxes=True,shared_yaxes=True)


        for count, gene in enumerate(genes):
            
            row = count // n_col
            col = count % n_col

            gene_tpm, gene_tpm_min, gene_tpm_max = input_arr[gene], input_arr[gene].min(), input_arr[gene].max()
            gene_tpm_norm = (gene_tpm - gene_tpm_min) / (gene_tpm_max - gene_tpm_min)

            scatter = go.Scatter(x = input_arr["x_umap1"], y = input_arr["x_umap2"],
                                    mode = 'markers',
                                    marker=dict(
                                        size=10,  # Adjust marker size if needed
                                        color=gene_tpm_norm,  # Color by the numerical value
                                        colorscale='Viridis',
                                        #colorbar=dict(title=gene,x=colorbar_x)#xanchor="right",yanchor="middle")  # Choose a colorscale (e.g., Viridis, Plasma, etc.)),
                                        ))
            
            fig.add_trace(scatter, row= row+1, col= col+1)
        
        fig.update_layout(height=300*n_rows, width=300*n_col, showlegend=False)
        fig.update_xaxes(scaleanchor="y", scaleratio=1, showticklabels=False)
        fig.update_yaxes(scaleanchor="x", scaleratio=1, showticklabels=False)

        return fig
    

def plot_pancan_archetypes():

    fig = px.scatter(sh.categorical_data, x = "x_umap1", y="x_umap2", 
                         color="archetype", color_discrete_map=sh.colors_pancan)
    fig.update_traces(marker=dict(size=12))
    fig.update_layout(legend_title_text = "Archetype")
    fig.update_layout(template = "simple_white",
                      legend=dict(
                        y=0.5,  # Center vertically
                        #xanchor='left',  # Align to the left of the legend box
                        #yanchor='middle',  # Center vertically
                        font=dict(size=18)  # Increase font size
                        ),
    )
    fig.update_yaxes(visible=False)
    fig.update_xaxes(visible=False)
    #fig.update_layout(width =1200)


    return fig