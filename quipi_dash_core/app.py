from shiny import App, render, ui
import plotly.express as px
from shinywidgets import output_widget, render_widget 

from shared import app_dir, df, genes, non_genes

import numpy as np
import pandas as pd
pd.options.plotting.backend = 'plotly'
import matplotlib.pyplot as plt

# Define the UI
app_ui = ui.page_fluid(
    ui.panel_title("QuIPI - Querying IPI"),
    
    ui.navset_tab(
        ui.nav_panel("Home", ui.h2("Welcome to the Home Page"), ui.p("This is the home tab.")),

        ui.nav_panel("PanCan UMAP Expression Overlay",
                  ui.input_selectize("pancan_gene_input",
                                     "Select Genes:",
                                     genes,
                                     multiple = True),
                  ui.output_ui("plots_multiple",
                               inline=True,
                               fill = True,
                               fillable=True)),
        ui.nav_panel("Barplots")
    ),

    title = "QuIPI - Querying IPI",
)

# Define the server logic (empty in this case since we're just creating static pages)
def server(input, output, session):
    
    def pancan_multiple_expression_plot(input_gene):
        @render_widget
        def pancan_single_gene_fig():

            input_arr = df[input_gene]
            logged = np.where(input_arr != 0, np.log2(input_arr), 0)

            fig = df.plot.scatter(x="X_umap1", y = "X_umap2", 
                                color = logged,
                                color_continuous_scale="Viridis")
            fig.update_layout(template='simple_white')
            fig.update_layout(title_text= input_gene + " Log(TPM)", title_x=0.5)
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

        print(plot_output)
        print(ui.TagList(plot_output))
        return ui.TagList(plot_output)


    '''
    def render_plot_func(j):
        @render.plot
        def f():
            fig = plt.plot(range(1, j + 1), range(1, j + 1))
            return fig
        return f
    
    @output
    @render.ui
    def plots():
        plot_output_list = []
        for i in range(1, input.n() + 1):
            plotname = f"plot{i}"

            plot_output_list.append(ui.output_plot(plotname))
            output(render_plot_func(i), id=plotname)
        return ui.TagList(plot_output_list)
    '''
    

# Create the Shiny app
app = App(app_ui, server)

# Run the app
#app.run()