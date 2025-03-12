from shiny import App, render, ui, reactive
from shinyswatch import theme
import plotly.express as px
from shinywidgets import output_widget, render_widget 

import shared as sh
import quipi_humu.flow_boxplot as bh
import gex_violin as gv

import numpy as np
import pandas as pd

from pathlib import Path

tabs_mapped_to_gene_inputs = {"Violin Plots" : ["gex_viol_gene"],
}


#RUN_STYLE="background-color: #AFE1AF; color: black;"

panel_color = "#f0f0f0"


# Define the UI
app_ui = ui.page_fluid(

    ui.tags.style("""
        body { background-color: #e5e7f8; }  /* Background color */
        .nav-link { font-size: 20px; }
    """),
    ui.head_content(
        ui.tags.link(
            rel="icon", type="image/png", sizes="64x64", href="favicon-32x32.png"
        ),
    ),
    ui.tags.div(
        ui.tags.div(
            ui.tags.img(src="humu.png", style="width: 100px; margin-left: 10px; padding-top: 10px"),  # Left-aligned image
            ui.tags.span("QuIPI - HuMu", style="font-size: 50px; font-weight: bold;"),
            style="display: flex; align-items: flex-end; gap: 10px;" # Flexbox for horizontal alignment
        ),
        style="padding-bottom: 5px;"
    ),

    ui.page_navbar(
        ui.nav_panel(
            "Home",
            ui.head_content(
            ),
                ui.panel_title("QuIPI - HuMu"),
                ui.p("Here are some useful references."),
        ),

        ui.nav_panel("Flow Boxplots",
            #ui.h4("Explore gene expression by category"),
            ui.layout_sidebar(
                ui.sidebar(
                    ui.input_selectize("box_score_1", "Choose First Score", sh.flow_scores),
                    ui.input_selectize("box_score_2", "Choose Second Score if ratio desired", ["---"] + sh.flow_scores, selected="---"),
                    ui.input_selectize("box_x_cat", "Select X-Axis Category", ["Species", "Group"]),
                    ui.input_selectize("box_x_cat_filter", "**Subset X-Axis Categories.**", [], multiple=True),
                    ui.input_selectize("box_group", "Group By:",  ["---"] + sh.flow_cats, selected="---"),
                    ui.input_action_button("box_run", "RUN"),
                    bg=panel_color
                ),
                ui.card(ui.card_body(output_widget("expression_box_viol")),
                        ui.card_footer("Click button in the bottom right for fullscreen view."),
                        full_screen=True),
                bg=panel_color,
            )
        
        ),

        ui.nav_panel("Violin Plots",
            ui.h4("Explore gene expression"),
            ui.layout_sidebar(
                ui.sidebar(
                    ui.input_selectize("gex_viol_gene", "Choose Gene to plot:", []),
                    ui.input_selectize("gex_viol_x_cat", "Choose X-Axis Category:", sh.categoricals_opts),
                    ui.input_selectize("gex_viol_cat_subset", "Subset Categories:", [], multiple=True),
                    ui.input_selectize("gex_viol_groupby", "Group by:", ["---"] + sh.categoricals_opts, selected="---"),
                    ui.input_selectize("gex_viol_splitby", "Split by:", ["---"] + sh.categoricals_opts, selected="---"),
                    ui.input_action_button("gene_viol_run", "RUN"),
                    bg=panel_color
                ),
                ui.card(output_widget("gex_viol"), full_screen=True),
                bg=panel_color
            )       
        ),
        ui.nav_spacer(),
        ui.nav_control(ui.a("Return to QuIPI", href="https://quipi.org/app/quipi", target="_blank", class_="nav-link")),
        id = "quipi_top_nav",
        theme=theme.cosmo,
        bg= '#1a1807',
    ),
)   

def server(input, output, session):

    ##### FLOW BOX PLOTS

    @reactive.effect
    @reactive.event(input.box_x_cat)  # Trigger when category changes
    def update_box_viol_selectize():
        x_cat = input.box_x_cat()
        flow_df = pd.read_feather("./quipi_humu_data/quipi_humu_flow_table.feather", columns=[x_cat])
        cats = list(flow_df[x_cat].unique())
        ui.update_selectize("box_x_cat_filter", choices=cats, selected=cats)
        
    @render_widget
    @reactive.event(input.box_run)
    def expression_box_viol():
        score1 = input.box_score_1()
        score2 = input.box_score_2()
        x_cat = input.box_x_cat()
        x_cat_filter = input.box_x_cat_filter()
        group = input.box_group()
        fig = bh.box_humu_flow(score1, score2, x_cat, x_cat_filter, group)
        return fig
    
    ##### GEX Violin Plots
    @render_widget
    @reactive.event(input.gene_viol_run)
    def gex_viol():
        gene = input.gex_viol_gene()
        x_cat = input.gex_viol_x_cat()
        x_sub = input.gex_viol_cat_subset()
        group = input.gex_viol_groupby()
        split = input.gex_viol_splitby()

        fig = gv.plot_sc_violin(gene,x_cat,x_sub,group,split)
        return fig

    @reactive.effect
    @reactive.event(input.gex_viol_x_cat)  # Trigger when category changes
    def update_box_viol_selectize():
        x_cat = input.gex_viol_x_cat()
        cat_opts = list(pd.read_feather("./quipi_humu_data/quipi_humu_adata_clean_full_PROC.feather", columns=[x_cat])[x_cat].unique())
        ui.update_selectize("gex_viol_cat_subset", choices=cat_opts, selected=cat_opts)
    


    ##### HELPER 

    tabs_visited = []
    ## HELPER: Populate individual gene selection boxes to avoid long startup.
    @reactive.effect
    def populate_gene_selections():
        current_tab = input.quipi_top_nav()  # Read the active tab
        if current_tab not in tabs_visited:
            tabs_visited.append(current_tab)
            if current_tab in tabs_mapped_to_gene_inputs:
                for id in tabs_mapped_to_gene_inputs[current_tab]:
                    ui.update_selectize(
                        id,
                        choices=sh.genes,
                        selected=[],
                        server=True,
                    )

# Create the Shiny app
app_dir = Path(__file__).parent
app = App(app_ui, server,static_assets= app_dir / "www")



