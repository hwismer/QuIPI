from shiny import App, render, ui, reactive
import plotly.express as px
from shinywidgets import output_widget, render_widget 
from plotly.subplots import make_subplots
import plotly.graph_objects as go

import shared as sh
#from shared import app_dir, quipi_raw, quipi_log10, quipi_log2
#from shared import genes, non_genes, categoricals, colors_pancan, colors_indic, transformations, corr_methods

import numpy as np
import pandas as pd
pd.options.plotting.backend = 'plotly'
import matplotlib.pyplot as plt


from scipy.stats import zscore
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from scipy.stats import ranksums


# Define the UI
app_ui = ui.page_fluid(
    ui.panel_title("QuIPI - Querying IPI"),
    
    ui.navset_tab(

        ui.nav_panel("Home", 
                     ui.h2("Welcome to the Home Page"), 
                     ui.p("This is the home tab.")),

        ui.nav_panel("PanCan UMAP Gene Expression",
                     
                    ui.page_sidebar(
                         ui.sidebar(
                    
                            ui.input_selectize("pancan_gene_input",
                                                "Select Genes:",
                                                sh.genes,
                                                multiple = True),
                            ui.input_selectize("pancan_compartment_input",
                                               "Select Compartment:",
                                               list(sh.quipi_raw["compartment"].unique())),
                            ui.input_selectize("pancan_umap_transformation",
                                            "Choose Transformation: ",
                                            ["Raw", "Log2", "Log10"],
                                            multiple= False,
                                            selected= "Log2"),
                            output_widget("pancan_archetypes")),

                        ui.layout_columns(output_widget("pancan_subplots")))),

        ui.nav_panel("Boxplots / Violin Plots",
                    ui.page_sidebar(
                        ui.sidebar(
                            ui.input_selectize("box_viol_plot",
                                               "Select Plot Type:",
                                               ["Boxplot","Violin Plot"]),
                            ui.input_selectize("box_viol_gene_input",
                                                "Select Gene:",
                                                sh.genes,),
                            ui.input_selectize("box_viol_x_category",
                                               "Select X-Axis Category:",
                                               list(sh.categoricals_dict.keys())),
                            ui.input_selectize("box_viol_groupby",
                                            "Group by:",
                                            list(sh.categoricals_dict.keys())),
                            ui.input_selectize("box_viol_transformation",
                                            "Choose Transformation: ",
                                            ["Raw", "Log2", "Log10"],
                                            multiple= False,
                                            selected= "Log2")),             
                        output_widget("expression_box_viol"))),

        ui.nav_panel("Correlation Matrix Analysis",
            ui.page_sidebar(
                ui.sidebar("Explore the correlation of selected genes.",
                ui.input_selectize("corr_gene_input",
                                    "Select Genes:",
                                    sh.genes,
                                    multiple=True),
                ui.input_selectize("corr_indication",
                                   "Select Indications:",
                                   list(sh.quipi_raw["indication"].unique()),
                                   multiple = True,
                                   selected=list(sh.quipi_raw["indication"].unique())),
                ui.input_selectize("corr_tissue",
                                   "Select Tissue:",
                                   ["Tumor", "Normal"],
                                   multiple = True,
                                   selected=["Tumor", "Normal"]),
                ui.input_selectize("corr_compartment",
                                   "Select Compartments:",
                                   list(sh.quipi_raw["compartment"].unique()),
                                   multiple=True,
                                   selected=list(sh.quipi_raw["compartment"].unique())),
                ui.input_selectize("corr_archetype",
                                   "Select Archetypes:",
                                   list(sh.quipi_raw["archetype"].unique()),
                                   multiple=True,
                                   selected=list(sh.quipi_raw["archetype"].unique())),
                ui.input_selectize("corr_transform",
                                   "Select Transformation:",
                                   ["Raw", "Log2", "Log10"],
                                   selected="Log2"),
                ui.input_selectize("corr_method_input",
                                "Select Method:",
                                ["Pearson", "Spearman"],
                                selected="Spearman")),
            output_widget("gene_correlation_heatmap"))),
            
            ui.nav_panel("Gene Factor Analysis",
                         ui.page_sidebar(
                             ui.sidebar(
                             ui.input_selectize("gene_factor_genes",
                                                "Select Genes:",
                                                sh.genes,
                                                multiple=True),
                             ui.input_selectize("gene_factor_compartment",
                                                "Select Compartment:",
                                                list(sh.quipi_raw["compartment"].unique())),
                             ui.input_slider("gene_factor_slider",
                                             "Choose Percentile:",
                                             value = .2,
                                             min = 0.05,
                                             max = .5,
                                             step = .05),
                                             ),
                             output_widget("gene_factor_analysis")))),
    title = "QuIPI - Querying IPI",
)

# Define the server logic (empty in this case since we're just creating static pages)
def server(input, output, session):

    @render_widget
    def pancan_subplots():
        
        transform = input.pancan_umap_transformation()
        genes = input.pancan_gene_input()
        compartment = input.pancan_compartment_input()
        input_arr = sh.transformations[transform]
        input_arr = input_arr[input_arr["compartment"] == compartment]

        if len(genes) != 0:

            n_col = min(4,len(genes))
            n_rows = (len(genes) + n_col - 1) // n_col

            fig = make_subplots(rows =  n_rows , cols = n_col, 
                                subplot_titles=genes,
                                vertical_spacing=.05,horizontal_spacing=.01,
                                shared_xaxes=True,shared_yaxes=True)


            for count, gene in enumerate(genes):
                row = count // n_col
                col = count % n_col

                scatter = go.Scatter(x = input_arr["x_umap1"], y = input_arr["x_umap2"],
                                        mode = 'markers',
                                        marker=dict(
                                            size=10,  # Adjust marker size if needed
                                            color=input_arr[gene],  # Color by the numerical value
                                            colorscale='Viridis',
                                            colorbar=dict(title=gene,xanchor="right",yanchor="middle")  # Choose a colorscale (e.g., Viridis, Plasma, etc.)),
                                            ))
                
                fig.add_trace(scatter, row= row+1, col= col+1)
            
            fig.update_layout(height=300*n_rows, width=300*n_col, showlegend=False)
            fig.update_xaxes(scaleanchor="y", scaleratio=1, showticklabels=False)
            fig.update_yaxes(scaleanchor="x", scaleratio=1, showticklabels=False)

            return fig
        
        
    def pancan_archetypes():

        
        colors = [sh.colors_pancan[classification] for classification in sh.pancan_only_raw["archetype"]]
        fig = go.Scatter(x = sh.pancan_only_raw["x_umap1"], y = sh.pancan_only_raw["x_umap2"],
                         mode = 'markers',
                         marker=dict(
                               size=10,
                               color=colors,
                              ),
                         hovertext=sh.pancan_only_raw["archetype"],
                         showlegend=True
                        )
        #fig = sh.quipi_raw.plot.scatter(x="x_umap1", y="x_umap2",
        #                            color = "archetype",
        #                            color_discrete_map = sh.colors_pancan)
        #fig.update_layout(template='simple_white')
        #fig.update_layout(title_text= "Archetype UMAP", title_x=0.5)
        #fig.update_layout(autosize=False, width=500, height=400)
        return fig

    @render_widget
    def expression_box_viol():

        transform = input.box_viol_transformation()
        x_cat = sh.categoricals_dict[input.box_viol_x_category()]

        input_arr = sh.transformations[transform]

        gene = input.box_viol_gene_input()
        group = sh.categoricals_dict[input.box_viol_groupby()]

        if input.box_viol_plot() == "Boxplot":
            fig = px.box(input_arr, x = x_cat, y = gene, color = group,
                         color_discrete_sequence=px.colors.qualitative.D3)
        else:
            fig = px.violin(input_arr, x = x_cat, y = gene, color = group,
                            color_discrete_sequence=px.colors.qualitative.D3)
            
        fig.update_layout(title_text= gene + " " + transform + "(TPM)", title_x=0.5)
        return fig
    

    @render.data_frame
    def gene_corr_df():
        genes = list(input.corr_gene_input())
        indications = input.corr_indication()
        method = sh.corr_methods[input.corr_method_input()]
        compartments = input.corr_compartment()
        archetypes = input.corr_archetype()
        tissues = [sh.tissue_dict[tis] for tis in input.corr_tissue()]

        transform = input.corr_transform()

        input_arr = sh.transformations[transform]
        input_arr = input_arr[input_arr["indication"].isin(indications)]
        input_arr = input_arr[input_arr["compartment"].isin(compartments)]
        input_arr = input_arr[input_arr["archetype"].isin(archetypes)]
        input_arr = input_arr[input_arr["sample_type_cat"].isin(tissues)]
        input_arr = input_arr[genes]

        return input_arr


    @render_widget
    def gene_correlation_heatmap():

        genes = list(input.corr_gene_input())
        indications = input.corr_indication()
        method = sh.corr_methods[input.corr_method_input()]
        compartments = input.corr_compartment()
        archetypes = input.corr_archetype()
        tissues = [sh.tissue_dict[tis] for tis in input.corr_tissue()]


        transform = input.corr_transform()

        input_arr = sh.transformations[transform]
        #input_arr = input_arr[genes]
        input_arr = input_arr[input_arr["indication"].isin(indications)]
        input_arr = input_arr[input_arr["compartment"].isin(compartments)]
        input_arr = input_arr[input_arr["archetype"].isin(archetypes)]
        input_arr = input_arr[input_arr["sample_type_cat"].isin(tissues)]
        input_arr = input_arr[genes]

        corr_df = input_arr.corr(method=method)

        # Create a mask for the upper triangle
        mask = np.triu(np.ones_like(corr_df, dtype=bool))

        # Apply the mask to the correlation matrix
        corr_lower_tri = corr_df.mask(mask)

        corr_lower_tri = corr_lower_tri.dropna(axis=0, how = "all")
        corr_lower_tri = corr_lower_tri.dropna(axis=1, how = "all")

        fig = px.imshow(corr_lower_tri.fillna(""),
                        zmin = -1, zmax =1,
                        color_continuous_scale = "RdBu_r",
                        text_auto=True)
        fig.update_layout(template='simple_white')
        fig.update_xaxes(showgrid=False, showline=False)
        fig.update_yaxes(showgrid=False, showline=False)

        return fig
    
    @render_widget
    def gene_factor_analysis():
        gene_set = list(input.gene_factor_genes())
        compartment = input.gene_factor_compartment()
        percentile = input.gene_factor_slider()

        log2_subset = sh.quipi_log2[sh.quipi_log2["archetype"] != "Unclassified"][sh.quipi_log2["compartment"] == compartment][sh.non_genes+gene_set]
        raw_subset = sh.quipi_raw[sh.quipi_raw["archetype"] != "Unclassified"][sh.quipi_raw["compartment"] == compartment][sh.non_genes+gene_set]

        z_subset = log2_subset[gene_set].apply(zscore)
        z_subset["factor_score"] = z_subset.mean(axis=1)
        log2_subset_full = log2_subset[sh.non_genes].merge(z_subset,left_index=True,right_index=True)

        fig = make_subplots(rows = 1 , cols = 2,
                                vertical_spacing=.05,horizontal_spacing=.01,
                                shared_xaxes=True,shared_yaxes=True,
                                subplot_titles=[str(gene_set), "PanCan Archetypes"])
        

        test = go.Scatter(x = log2_subset_full["x_umap1"], y = log2_subset_full["x_umap2"],
                          mode = 'markers',
                          marker=dict(
                              size=10,
                              color=log2_subset_full["factor_score"],
                              colorscale='Viridis', 
                              ))
        #gene_factors.update_layout(template='simple_white')
        #gene_factors.update_layout(title_text= str(gene_set), title_x=0.5)

        fig.add_trace(test, row=1,col=1)
        fig.add_trace(pancan_archetypes(), row=1,col=2)
        fig.update_xaxes(scaleratio=1, showticklabels=False)
        fig.update_yaxes(scaleratio=1, showticklabels=False)
        fig.update_layout(showlegend=False)

        return fig


        '''
        n_total = len(log2_subset_full)
        num_tailed = int(n_total * percentile)

        log2_sorted = log2_subset_full.sort_values("factor_score",ascending=False)
        top = log2_sorted.head(num_tailed)
        bot = log2_sorted.tail(num_tailed)

        top_data = top[gene_set]
        bot_data = bot[gene_set]

        p_vals = ranksums(top_data,bot_data)[1]
        p_adj = multipletests(p_vals, method = "fdr_bh")[1]

        raw_top = raw_subset.loc[top.index][gene_set]
        raw_bot = raw_subset.loc[bot.index][gene_set]

        top_avg_tpm = raw_top.mean(axis=0) + .01
        bot_avg_tpm = raw_bot.mean(axis=0) + .01

        fc = np.log2(top_avg_tpm / bot_avg_tpm)

        fig = px.scatter(x = fc, 
           y = -np.log10(p_adj),
           text=gene_set)

        fig.update_layout(autosize=False, width=600, height=600)
        fc_abs_max = max(fc,key=abs)
        fig.update_xaxes(range=[-fc_abs_max - 1, fc_abs_max + 1])
        fig.update_layout(template='simple_white')
        '''

# Create the Shiny app
app = App(app_ui, server)

# Run the app
#app.run()


