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

TAG_STYLE="""
            /* Increase font size for nav panel text and center it */
            .nav-link {
                font-size: 20px;
                font-weight: bold;
                text-align: center;
                display: flex;
                align-items: center;
                justify-content: center;
            }

            /* Ensure active tab text is also centered */
            .nav-link.active {
                text-align: center;
                display: flex;
                align-items: center;
                justify-content: center;
            }

            /* Optional: Add padding for better spacing */
            .nav-link {
                padding: 20px;
            }
            """

# Define the UI
app_ui = ui.page_navbar(

    ui.nav_panel(
        "Home",
        ui.head_content(
            ui.tags.style(
                TAG_STYLE
            )
        ),
            ui.panel_title("QuIPI - HuMu"),
            ui.p("Here are some useful references."),
    ),

    ui.nav_panel("Flow Boxplots",
        ui.h4("Explore gene expression by category"),
        ui.layout_sidebar(
            ui.sidebar(
                ui.input_selectize("box_score_1", "Choose First Score", sh.flow_scores),
                ui.input_selectize("box_score_2", "Choose Second Score if ratio desired", ["---"] + sh.flow_scores, selected="---"),
                ui.input_selectize("box_x_cat", "Select X-Axis Category", ["Species", "Group"]),
                ui.input_selectize("box_x_cat_filter", "**Subset X-Axis Categories.**", [], multiple=True),
                ui.input_selectize("box_group", "Group By:",  ["---"] + sh.flow_cats, selected="---"),
                ui.input_action_button("box_run", "RUN")     
            ),
        ui.card(ui.card_body(output_widget("expression_box_viol")),
                ui.card_footer("Click button in the bottom right for fullscreen view."),
                full_screen=True)
        )
    ),

    ui.nav_spacer(),
    ui.nav_control(
        ui.tags.a(
        ui.tags.i(class_="fa fa-github", style="font-size: 40px; color: black; justify-content: center; margin-top: 8px"),
        href="https://github.com/HarrisonWismer/QuIPI",  # The URL to navigate to
        target="_blank",  # Opens in a new tab
        style="text-decoration: none; margin-left: 0px; justify-content: center;",  # Optional styling
        ),
    ),
    id = "quipi_top_nav",
    title = ui.div(ui.tags.a(ui.img(src="humu.png", height="65px"))),
    theme=theme.lumen,
    bg= '#85aad4',
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



