from shiny import App, render, ui, reactive
from shiny.types import ImgData
from shinyswatch import theme
import plotly.express as px
from shinywidgets import output_widget, render_widget 

import shared as sh
import box_viol_expression_humu as bh

import numpy as np
import pandas as pd
from io import StringIO

from pathlib import Path

import asyncio
import scipy

RUN_STYLE="background-color: #AFE1AF; color: black;"

panel_color = "#f0f0f0"


# Define the UI
app_ui = ui.page_fluid(
    ui.tags.style("""
        body { background-color: #e5e7f8; }  /* Background color */
        .nav-link { font-size: 20px; }
    """),

    ui.tags.div(
        ui.tags.div(
            ui.tags.img(src="humu_edit2.png", style="width: 100px; margin-left: 10px; padding-top: 10px"),  # Left-aligned image
            ui.tags.span("QuIPI - HuMu", style="font-size: 50px; font-weight: bold;"),
            style="display: flex; align-items: flex-end; gap: 10px;" # Flexbox for horizontal alignment
        ),
        style="padding-bottom: 5px;"
    ),

    #ui.a("Return to QuIPI", href="https://quipi.org/app/quipi", target="_blank", class_="nav-link"),

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
                    ui.input_selectize("viol_genes", "Choose Genes to plot:", sh.adata_vars),
                    #ui.input_selectize("")
                    bg=panel_color
                ),
                bg=panel_color

            )       
        ),
        ui.nav_spacer(),
        #ui.input_action_button("go_to_site", "Visit Website", class_="btn-primary"),
        #ui.nav_panel("Return to QuIPI"),
        ui.nav_control(ui.a("Return to QuIPI", href="https://quipi.org/app/quipi", target="_blank", class_="nav-link")),
        id = "quipi_top_nav",
        theme=theme.cosmo,
        bg= '#1a1807',
    )
)   

def server(input, output, session):

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

        fig = bh.box_humu(score1, score2, x_cat, x_cat_filter, group)

        return fig

# Create the Shiny app
app_dir = Path(__file__).parent
app = App(app_ui, server,static_assets= app_dir / "www")



