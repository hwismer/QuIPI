from shiny import App, render, ui, reactive
from shiny.types import ImgData
from shinyswatch import theme
import plotly.express as px
from shinywidgets import output_widget, render_widget 

#import shared as sh

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
            ui.panel_title("Welcome to QuIPI - HuMu"),
            ui.p("Here are some useful references."),
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
    pass

# Create the Shiny app
app_dir = Path(__file__).parent
app = App(app_ui, server,static_assets= app_dir / "www")



