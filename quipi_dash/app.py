#import seaborn as sns

import pandas as pd
import numpy  as np

import plotly.io as pio
pio.renderers.default = 'iframe'
import plotly.express as px
pd.options.plotting.backend = 'plotly'

#from faicons import icon_svg

# Import data from shared.py
from shared import app_dir, df, genes, non_genes

from shiny import reactive
from shiny.express import input, render, ui
from shinywidgets import render_plotly

ui.page_opts(title="QuIPI - Querying IPI", fillable=True)

with ui.navset_tab(id = "primary_tab"):
    with ui.nav_panel("PanCan Archetype Feature Plot"):
        "Display Features on PanCan UMAP"

        ui.input_selectize(  
        "pancan_umap_overlay_genes",  
        "Select options below:",
        genes,
         multiple=True,)  
        
        @render.text
        def value():
            return f"{input.selectize()}"


        with ui.card(full_screen=True):
            ui.card_header("PanCan UMAP")

            # np.where(df != 0, np.log2(df),0)
            @render_plotly
            def pancan_expression_plot():

                input_gene = input.pancan_umap_overlay_genes()[0]
                input_arr = df[input_gene]
                logged = np.where(input_gene != 0, np.log2(input_arr), 0)

                fig = df.plot.scatter(x="X_umap1", y = "X_umap2", 
                                    color = logged,
                                    color_continuous_scale="Viridis")
                fig.update_layout(template='simple_white')
                fig.update_layout(title_text= input_gene + " Log(TPM)", title_x=0.5)
                fig.update_layout(autosize=False, width=500, height=400)
                return fig

    with ui.nav_panel("Scatter"):
        "Panel B Content"

    with ui.nav_panel("Heatmap"):
        "Panel B Content"

    with ui.nav_panel("Boxplot"):
        "Panel B Content"

    with ui.nav_panel("Violin Plot"):
        "Panel B Content"

    with ui.nav_panel("Barplot"):
        "Panel B Content"
    
    

#with ui.sidebar(title="Select Genes"):
#    ui.input_slider("mass", "Mass", 2000, 6000, 6000)
#    ui.input_checkbox_group(
#        "species",
#        "Species",
#        ["Adelie", "Gentoo", "Chinstrap"],
#        selected=["Adelie", "Gentoo", "Chinstrap"],
#    )

ui.include_css(app_dir / "styles.css")


@reactive.calc
def filtered_df():
    filt_df = df[df["species"].isin(input.species())]
    filt_df = filt_df.loc[filt_df["body_mass_g"] < input.mass()]
    return filt_df

