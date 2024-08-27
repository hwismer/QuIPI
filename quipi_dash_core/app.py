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


# Define the UI
app_ui = ui.page_fluid(
    ui.panel_title("QuIPI - Querying IPI"),

    #ui.layout_columns(output_widget("pancan_archetypes"),output_widget("pancan_archetypes2"),
    #                      width="300px", gap="20px",justify="center"), 
    
    ui.navset_tab(

        ui.nav_panel("Home", ui.h2("Welcome to the Home Page"), ui.p("This is the home tab.")),

        ui.nav_panel("PanCan UMAP Expression Overlay",
                     
                     ui.page_sidebar(
                         ui.sidebar(
                    
                            ui.input_selectize("pancan_gene_input",
                                                "Select Genes:",
                                                sh.genes,
                                                multiple = True),
                            ui.input_selectize("pancan_umap_transformation",
                                            "Choose Transformation: ",
                                            ["Raw", "Log2", "Log10"],
                                            multiple= False,
                                            selected= "Log2")),

                        ui.layout_columns(ui.output_ui("plots_multiple"),
                                                width = "300px", gap = "20px", justify="center"))),
        ui.nav_panel("PanCan UMAP Testing",
                     
                    ui.page_sidebar(
                         ui.sidebar(
                    
                            ui.input_selectize("pancan_gene_input_2",
                                                "Select Genes:",
                                                sh.genes,
                                                multiple = True),
                            ui.input_selectize("pancan_umap_transformation_2",
                                            "Choose Transformation: ",
                                            ["Raw", "Log2", "Log10"],
                                            multiple= False,
                                            selected= "Log2")),

                        ui.layout_columns(output_widget("pancan_subplots")))),

                  #output_widget("pancan_archetypes"),
                  #output_widget("pancan_archetype2"),
                  #ui.output_ui("plots_multiple",
                  #             inline=True,
                  #             fill = True,
                  #             fillable=True)),

        ui.nav_panel("Boxplots / Violin Plots",
                    ui.page_sidebar(
                        ui.sidebar(
                            ui.input_selectize("boxplot_gene_input",
                                                "Select Gene:",
                                                sh.genes,),
                            ui.input_selectize("boxplot_x_category",
                                               "Select X-Axis Category:",
                                               list(sh.categoricals_dict.keys())),
                            ui.input_selectize("boxplot_groupby",
                                            "Group by:",
                                            list(sh.categoricals_dict.keys())),
                            ui.input_selectize("boxplot_transformation",
                                            "Choose Transformation: ",
                                            ["Raw", "Log2", "Log10"],
                                            multiple= False,
                                            selected= "Log2")),             
                        output_widget("expression_boxplot"))),

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
            output_widget("gene_correlation_heatmap")))),
            #ui.output_data_frame("gene_corr_df")))),
    title = "QuIPI - Querying IPI",
)

# Define the server logic (empty in this case since we're just creating static pages)
def server(input, output, session):

    @render_widget
    def pancan_subplots():
        
        transform = input.pancan_umap_transformation_2()
        genes = input.pancan_gene_input_2()
        input_arr = sh.transformations[transform]

        n_col = 4

        n_rows = len(genes) % n_col

        single_subplot = 400


        fig = make_subplots(rows =  n_rows+1 , cols = n_col, subplot_titles=genes,
                            column_widths=[.25,.25,.25,.25],
                            row_heights=[1/(n_rows+1) for x in range(n_rows+1)],
                            shared_xaxes=True,
                            shared_yaxes=True)
                            #shared_xaxes=True, shared_yaxes=True,)
                            #column_widths=[1/n_col for x in range(n_col)],
                            #row_heights = [1/(n_rows+1) for x in range(n_rows+1)])

        for count, gene in enumerate(genes):
            row = count // n_col
            col = count % n_col

            scatter = go.Scatter(x = input_arr["x_umap1"], y = input_arr["x_umap2"],
                                     mode = 'markers',
                                     marker=dict(
                                         size=10,  # Adjust marker size if needed
                                         color=input_arr[gene],  # Color by the numerical value
                                         colorscale='Viridis',  # Choose a colorscale (e.g., Viridis, Plasma, etc.)
                                         showscale=True),
                                         )
            
            fig.add_trace(scatter, row= row+1, col= col+1)

        #fig.update_layout(autosize=False, width = 1200, height = 300 * (n_rows+1),
        #                  margin=dict(l=00, r=00, t=00, b=00))

        return fig
    
    
    @render_widget
    def pancan_archetypes():
        fig = sh.quipi_raw.plot.scatter(x="x_umap1", y="x_umap2",
                                     color = "archetype",
                                     color_discrete_map = colors_pancan)
        fig.update_layout(template='simple_white')
        fig.update_layout(title_text= "Archetype UMAP", title_x=0.5)
        fig.update_layout(autosize=False, width=500, height=400)
        return fig
    
    @render_widget
    def pancan_archetypes2():
        fig = sh.quipi_raw.plot.scatter(x="x_umap1", y="x_umap2",
                                     color = "archetype",
                                     color_discrete_map = colors_pancan)
        #fig.update_layout(template='simple_white')
        fig.update_layout(title_text= "Archetype UMAP", title_x=0.5)
        fig.update_layout(autosize=False, width=500, height=400)
        return fig

    
    def pancan_multiple_expression_plot(input_gene):
        @render_widget
        def pancan_single_gene_fig():
            
            transform = input.pancan_umap_transformation()
            input_arr = sh.transformations[transform]

            #logged = np.where(input_arr != 0, np.log2(input_arr), 0)

            fig = input_arr.plot.scatter(x="x_umap1", y = "x_umap2", 
                                color = input_gene,
                                color_continuous_scale="Viridis",)
                                #hover_data = "archetype")
            fig.update_layout(template='simple_white')
            fig.update_layout(title_text= input_gene + " " + transform + "(TPM)", title_x=0.5)
            fig.update_layout(autosize=False, width=500, height=400)
            return fig
            
    
        return pancan_single_gene_fig
    
    
    
    @render.ui
    def plots_multiple():
        plot_output = []
        for gene in input.pancan_gene_input():
            plotname = f"plot{gene.replace(".","_")}"
            plot_output.append(output_widget(plotname))
            output(pancan_multiple_expression_plot(gene), id = plotname)

        return ui.TagList(plot_output)
    



    @render_widget
    def expression_boxplot():

        transform = input.boxplot_transformation()
        x_cat = sh.categoricals_dict[input.boxplot_x_category()]

        input_arr = sh.transformations[transform]

        gene = input.boxplot_gene_input()
        group = sh.categoricals_dict[input.boxplot_groupby()]

        fig = px.box(input_arr, x = x_cat, y = gene, color = group,
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






    

# Create the Shiny app
app = App(app_ui, server)

# Run the app
#app.run()


