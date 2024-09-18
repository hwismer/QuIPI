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

import molly_dge_split as mds


# Define the UI
app_ui = ui.page_fluid(
    ui.panel_title("QuIPI - Querying IPI"),
    
    ui.navset_tab(

        ui.nav_panel("Home", 
                     ui.h2("Welcome to QuIPI"), 
                     ui.p("Here are some useful references."),
                     ui.layout_column_wrap(
                        output_widget("pancan_archetypes"),
                        output_widget("cancer_glossary"))),

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
                                            selected= "Log2")),
                        ui.layout_columns(output_widget("pancan_subplots")))),

        ui.nav_panel("Box/Violin Plots",
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

        ui.nav_panel("Correlation Plots",
            ui.page_sidebar(
                    ui.sidebar("Explore the correlation of selected genes.",
                    ui.input_action_button("corr_run", "Run"),
                    ui.input_selectize("corr_gene_input",
                                        "Select Genes:",
                                        sh.genes,
                                        multiple=True),
                    ui.input_selectize("corr_indication",
                                    "Select Indications:",
                                    list(sh.quipi_raw["indication"].unique()),
                                    multiple = True),
                    ui.input_selectize("corr_tissue",
                                    "Select Tissue:",
                                    ["Tumor", "Normal"],
                                    multiple = True),
                    ui.input_selectize("corr_compartment",
                                    "Select Compartments:",
                                    list(sh.quipi_raw["compartment"].unique()),
                                    multiple=True),
                    ui.input_selectize("corr_archetype",
                                    "Select Archetypes:",
                                    list(sh.quipi_raw["archetype"].unique()),
                                    multiple=True),
                    ui.input_selectize("corr_transform",
                                    "Select Transformation:",
                                    ["Raw", "Log2", "Log10"],
                                    selected="Log2"),
                    ui.input_selectize("corr_method_input",
                                    "Select Method:",
                                    ["Pearson", "Spearman"],
                                    selected="Spearman")),
            output_widget("gene_correlation_heatmap"))),
            
            ui.nav_panel("PanCan Gene Factor Scores",
                         ui.page_sidebar(
                             ui.sidebar(
                                ui.input_action_button("gene_factor_run", "Run"),
                             
                                ui.input_selectize("gene_factor_genes",
                                                    "Select Genes:",
                                                    sh.genes,
                                                    multiple=True),
                                ui.input_selectize("gene_factor_compartment",
                                                    "Select Compartment:",
                                                    list(sh.quipi_raw["compartment"].unique()))),
                             output_widget("gene_factor_analysis"))),

            ui.nav_panel("Molly DGE Ranking",
                         ui.page_sidebar(
                             ui.sidebar(
                                ui.input_selectize("flow_score_to_rank",
                                                    "Feature Score For Ranking:",
                                                    list(sh.flow_scores.keys())),
                                ui.input_slider("dge_slider",
                                                "Choose Quartile:",
                                                value = .2,
                                                min = 0.05,
                                                max = .5,
                                                step = .05),
                                ui.input_selectize("dge_compartment",
                                                    "DGE Compartment:",
                                                    list(sh.quipi_raw["compartment"].unique())),
                             ),
                             ui.layout_column_wrap(
                                output_widget("compartment_featurescore_dge"),
                                #ui.output_table("compart_feat_table"),
                                #ui.output_table("compart_feat_table_bot")
    )))),
    title = "QuIPI - Querying IPI"
)

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
                                vertical_spacing=.05,horizontal_spacing=.02,
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
                                            #colorbar=dict(title=gene,x=(1+col)*.4756,y = (1+row) * .5)#xanchor="right",yanchor="middle")  # Choose a colorscale (e.g., Viridis, Plasma, etc.)),
                                            ))
                
                fig.add_trace(scatter, row= row+1, col= col+1)
            
            fig.update_layout(height=300*n_rows, width=300*n_col, showlegend=False)
            fig.update_xaxes(scaleanchor="y", scaleratio=1, showticklabels=False)
            fig.update_yaxes(scaleanchor="x", scaleratio=1, showticklabels=False)

            return fig
        
    @render_widget
    def cancer_glossary():
        df = sh.cancer_glossary_df

        fig = go.Figure(data=[go.Table(
            header=dict(values=list(df.columns),
                        fill_color='white',
                        font = dict(color = "black",size = 18),
                        align='center'),
            cells=dict(values=[df[col] for col in df.columns],
                    fill_color=[[sh.colors_indic[color] for color in df["Abbreviation"]]],  # Apply row colors
                    align='center',
                    height=30,
                    font = dict(color = 'white', size = 18)))
        ])
        
        fig.update_layout(autosize=False, width=600,height=800)
        return fig
        
    @render_widget
    def pancan_archetypes():
        

        fig = px.scatter(sh.pancan_only_raw, x = "x_umap1", y="x_umap2", 
                         color="archetype", color_discrete_map=sh.colors_pancan)
        fig.update_traces(marker=dict(size=12))
        fig.update_layout(legend_title_text = "Archetype")
        fig.update_layout(autosize=False, width=650, height=500,template = "simple_white")
        fig.update_yaxes(visible=False)
        fig.update_xaxes(visible=False)
   
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
                         color_discrete_sequence=px.colors.qualitative.D3,
                         labels=sh.categoricals_dict_reversed)
        else:
            fig = px.violin(input_arr, x = x_cat, y = gene, color = group,
                            color_discrete_sequence=px.colors.qualitative.D3,
                            labels=sh.categoricals_dict_reversed)
            
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
    @reactive.event(input.corr_run)
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
        fig.update_layout(template='simple_white', coloraxis_colorbar_x=.8)
        fig.update_xaxes(showgrid=False, showline=False)
        fig.update_yaxes(showgrid=False, showline=False)

        return fig
    
    @render_widget
    @reactive.event(input.gene_factor_run)
    def gene_factor_analysis():
        gene_set = list(input.gene_factor_genes())
        compartment = input.gene_factor_compartment()
        #percentile = input.gene_factor_slider()

        log2_subset = sh.quipi_log2[sh.quipi_log2["archetype"] != "Unclassified"][sh.quipi_log2["compartment"] == compartment][sh.non_genes+gene_set]
        raw_subset = sh.quipi_raw[sh.quipi_raw["archetype"] != "Unclassified"][sh.quipi_raw["compartment"] == compartment][sh.non_genes+gene_set]

        z_subset = log2_subset[gene_set].apply(zscore)
        z_subset["factor_score"] = z_subset.mean(axis=1)
        log2_subset_full = log2_subset[sh.non_genes].merge(z_subset,left_index=True,right_index=True)

        fig = px.scatter(log2_subset_full, x = "x_umap1", y = "x_umap2", 
                         color = "factor_score", color_continuous_scale="viridis")

        fig.update_layout(template="simple_white", autosize=False, width=650, height=500,
                          legend_title_text = "Factor Score")
        return fig

    
    @render_widget
    def compartment_featurescore_dge():

        return mds.molly_dge(input.flow_score_to_rank(), input.dge_compartment(), input.dge_slider())[0]
    
    @render.table
    def compart_feat_table():
        df,top,bot = mds.molly_dge(input.flow_score_to_rank(), input.dge_compartment(), input.dge_slider())
        #print(len(df),df["patient"], type(df["patient"]))
        #return pd.DataFrame.from_dict({"Test":[1,2]})
        #return mds.molly_dge(input.flow_score_to_rank(), input.dge_compartment(), input.dge_slider())[1]
        return top
    
    @render.table
    def compart_feat_table_bot():
        df,top,bot = mds.molly_dge(input.flow_score_to_rank(), input.dge_compartment(), input.dge_slider())
        #print(len(df),df["patient"], type(df["patient"]))
        #return pd.DataFrame.from_dict({"Test":[1,2]})
        #return mds.molly_dge(input.flow_score_to_rank(), input.dge_compartment(), input.dge_slider())[1]
        return bot
        '''
        rank_cat = sh.flow_scores[input.flow_score_to_rank()]
        comp = input.dge_compartment()
        quantile = input.dge_slider()

        top_patients = sh.quipi_flow[sh.quipi_flow[rank_cat] > sh.quipi_flow[rank_cat].quantile(1 - quantile)]["sample_name"]
        bot_patients = sh.quipi_flow[sh.quipi_flow[rank_cat] < sh.quipi_flow[rank_cat].quantile(quantile)]["sample_name"]


        top_tpm = sh.quipi_raw[sh.quipi_raw["sample_name"].isin(top_patients)]
        top_comp = top_tpm[top_tpm["compartment"] == comp]

        bot_tpm = sh.quipi_raw[sh.quipi_raw["sample_name"].isin(bot_patients)]
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

        fig.update_layout(autosize=False, width=600, height=600)
        fc_abs_max = abs(max(fc,key=abs))
        fig.update_layout(xaxis=dict(range=[-fc_abs_max-1, fc_abs_max+1]))
        fig.update_layout(template='simple_white')

        return fig

        

        
        ##########
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


