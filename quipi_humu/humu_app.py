from shiny import App, render, ui, reactive
from shinyswatch import theme
import plotly.express as px

px.defaults.template = "simple_white"
from shinywidgets import output_widget, render_widget 

#import quipi_shared as qsh
import humu_shared as hsh
import pandas as pd

from pathlib import Path

import humu_flow_boxplot as hfb
import humu_gex_plots as hxp

humu_tabs_mapped_to_gene_inputs = {"Boxplots" : [["humu_gex_box_gene", "Mouse"]],
                                   "Dotplots" : [["humu_gex_dot_gene", "Mouse"]],
                                   "HuMu Expression Comparison" : [["humu_box_comp_human_genes", "Human"],["humu_box_comp_mu_genes", "Mouse"]],
}


#RUN_STYLE="background-color: #AFE1AF; color: black;"

panel_color = "#f0f0f0"


# Define the UI
app_ui = ui.page_fluid(

    ui.tags.style("""
        body { background-color: #e5e7f8; }  /* Background color */
        .nav-link { font-size: 20px;
                    color: white;
        }
        .navbar-nav .nav-link.active {
            color: white !important;  /* Keep text white when selected */
        }
    """),
    ui.head_content(
        ui.tags.link(
            rel="icon", type="image/png", sizes="64x64", href="favicon-32x32.png"
        ),
    ),
    ui.tags.div(
        ui.tags.div(
            ui.tags.img(src="humu.png", style="width: 90px; margin-left: 10px; margin-top: 18px; margin-bottom: 8px",class_="my-image"),  # Left-aligned image
            ui.tags.span("QuIPI" + "\n" + "HuMu", style="font-size: 35px; font-weight: bold;"),
            style="display: flex; align-items: flex-end; gap: 10px;" # Flexbox for horizontal alignment
        ),
        style="padding-bottom: 5px;"
    ),
    ui.page_navbar(

        ui.nav_panel("Home",
            ui.layout_column_wrap(
                ui.card(
                    ui.h3("The Human-to-Mouse Cancer Translator Project (HuMu) aims to immuno-profile a series of \
                                common and exceptional models of cancer in mice to benchmark them against the diversity of \
                                TMEs in Human cancer (i.e immune archetypes described in Combes, Samad, et al. Cell 2022).")
                ),

                ui.card(ui.tags.img(src="quipi_humu_reference.png", width="450px")),

                ui.card(ui.h3("The HuMu dataset contains high-dimensional cytometry data (CyTOF) of 15 murine models that \
                    we use to study high level tumor-immune composition, as well as single-cell sequencing data \
                    from 9 of these models that we use to dissect more granular gene expression profiles across \
                    populations and tumor models."), style="display: flex; align-items: center; justify-content: center;"
                )
            )
        ),
        ui.nav_panel("HuMu Expression Comparison",
            ui.layout_sidebar(
                ui.sidebar(
                    ui.card(
                        ui.input_selectize("humu_box_comp_human_genes", "Choose Human Gene", [], multiple=False),
                        ui.output_ui("humu_box_viol_multiple_genes"),
                        ui.input_selectize("humu_box_comp_x_cat_filter",
                                            "**Subset Compartments**",
                                            hsh.quipi_compartments,
                                            selected=hsh.humu_compartments,
                                            multiple=True,
                                            remove_button=True,options={"plugins": ["clear_button"]}),
                        ui.input_selectize("humu_box_comp_transformation",
                                            "Select TPM Transformation: ",
                                            ["TPM", "Log2(TPM)"],
                                            multiple= False,
                                            selected= "Log2(TPM)"),
                    ),
                    ui.card(
                        ui.input_selectize("humu_box_comp_mu_genes", "Choose Murine Gene", []),
                        ui.input_selectize("humu_box_comp_cat_mouse_subset", "Subset Categories:",hsh.humu_compartments, selected=hsh.humu_compartments, multiple=True, remove_button=True,options={"plugins": ["clear_button"]}),
                    ),
                    ui.input_action_button("humu_box_comp_run", "RUN"),
                bg=panel_color
                ),
                ui.card(
                    ui.card(ui.card_header("Human"),
                            output_widget("humu_gene_comparison_human")),
                    ui.card(ui.card_header("Mouse"),
                            output_widget("humu_gene_comparison_mouse")),
                    full_screen=True
                ),
                bg=panel_color
            )
        ),

        ui.nav_menu("Mouse Gene Expression",

            ui.nav_panel("Boxplots",
                ui.layout_sidebar(
                    ui.sidebar(
                        ui.input_selectize("humu_gex_box_gene", "Choose Gene to plot:", []),
                        ui.input_selectize("humu_gex_box_x_cat", "Choose X-Axis Category:", hsh.categoricals_opts),
                        ui.input_selectize("humu_gex_box_cat_subset", "Subset Categories:", [], multiple=True, remove_button=True,options={"plugins": ["clear_button"]}),
                        ui.input_selectize("humu_gex_box_groupby", "Group by:", ["---"] + hsh.categoricals_opts, selected="---"),
                        ui.input_selectize("humu_gex_box_splitby", "Split by:", ["---"] + hsh.categoricals_opts, selected="---"),
                        ui.input_action_button("humu_gene_box_run", "RUN"),
                        bg=panel_color
                    ),
                    ui.card(output_widget("gex_box"), full_screen=True),
                    bg=panel_color
                )       
            ),

            ui.nav_panel("Dotplots",
                ui.layout_sidebar(
                    ui.sidebar(
                        ui.input_selectize("humu_gex_dot_gene", "Choose Genes to plot:", [], multiple=True),
                        ui.input_selectize("humu_gex_dot_groupby", "Group By:", hsh.categoricals_opts, selected="---"),
                        ui.input_selectize("humu_gex_dot_groups", "Subset Groupby Categories:", [], multiple=True),
                        ui.input_selectize("humu_gex_dot_splitby", "Split By:", ["---"] + hsh.categoricals_opts, selected="---"),
                        ui.input_selectize('humu_gex_dot_splits', "Subset Splitby Categories:", [], multiple=True),
                        ui.input_switch("humu_gex_dot_swap", "Swap Axes"),
                        ui.input_action_button("humu_gex_dot_run", "RUN"),
                        bg=panel_color
                    ),
                ui.card(ui.output_plot("humu_plot_gex_dotplot"), full_screen=True)
                )       
            ),
        ),

        ui.nav_panel("Flow-Score Boxplots",
            #ui.h4("Explore gene expression by category"),
            ui.layout_sidebar(
                ui.sidebar(
                    ui.input_selectize("humu_box_score_1", "Choose Parameter", hsh.flow_scores),
                    ui.input_selectize("humu_box_score_2", "Divide By Second Parameter (Optional)", ["---"] + hsh.flow_scores, selected="---"),
                    ui.input_switch("humu_box_switch", "Sort X-Category On Value"),
                    ui.input_selectize("humu_box_x_cat", "Select X-Axis Category", hsh.flow_cats),
                    ui.input_selectize("humu_box_x_cat_filter", "**Subset X-Axis Categories.**", [], multiple=True),
                    ui.input_selectize("humu_box_group", "Group By:",  ["---"] + hsh.flow_cats, selected="---"),
                    ui.input_action_button("humu_box_run", "RUN"),
                    bg=panel_color
                ),
                ui.card(ui.card_body(output_widget("humu_expression_box_viol")),
                        ui.card_footer("Click button in the bottom right for fullscreen view."),
                        full_screen=True),
                bg=panel_color,
            )
        ),
    ui.nav_spacer(),
    #ui.nav_control(ui.a("QuIPI HuMu", href="https://quipi.org/app/quipi", class_="nav-link")),
    id = "humu_top_nav",
    theme=theme.cosmo,
    bg = "#1a1807"
    )
)

def server(input, output, session):

    humu_tabs_visited = []

    ##### HUMU

    @reactive.effect
    @reactive.event(input.humu_box_x_cat)  # Trigger when category changes
    def update_box_viol_selectize():
        x_cat = input.humu_box_x_cat()
        flow_df = pd.read_feather("./quipi_humu_data/quipi_humu_flow_table.feather", columns=[x_cat])
        cats = list(flow_df[x_cat].unique())
        ui.update_selectize("humu_box_x_cat_filter", choices=cats, selected=cats)
        
    @render_widget
    @reactive.event(input.humu_box_run)
    def humu_expression_box_viol():
        score1 = input.humu_box_score_1()
        score2 = input.humu_box_score_2()
        x_cat = input.humu_box_x_cat()
        x_cat_filter = input.humu_box_x_cat_filter()
        group = input.humu_box_group()
        sort = input.humu_box_switch()
        fig = hfb.box_humu_flow(score1, score2, x_cat, x_cat_filter, group, sort)
        return fig
    
    ##### GEX Violin Plots
    @render_widget
    @reactive.event(input.humu_gene_box_run)
    def gex_box():
        gene = input.humu_gex_box_gene()
        x_cat = input.humu_gex_box_x_cat()
        x_sub = input.humu_gex_box_cat_subset()
        group = input.humu_gex_box_groupby()
        split = input.humu_gex_box_splitby()

        fig = hxp.plot_sc_box(gene,x_cat,x_sub,group,split)
        return fig

    @reactive.effect
    @reactive.event(input.humu_gex_box_x_cat)  # Trigger when category changes
    def humu_update_box_viol_selectize():
        x_cat = input.humu_gex_box_x_cat()
        cat_opts = list(pd.read_feather("./quipi_humu_data/quipi_humu_adata_clean_full_PROC.feather", columns=[x_cat])[x_cat].unique())
        ui.update_selectize("humu_gex_box_cat_subset", choices=cat_opts, selected=cat_opts)
    
    ##### GEX DOTPLOTS
    @reactive.effect
    @reactive.event(input.humu_gex_dot_groupby)
    def humu_update_dotplot_groupby():
        x_cat = input.humu_gex_dot_groupby()
        cat_opts = hsh.categorial_opts_dict[x_cat]
        ui.update_selectize("humu_gex_dot_groups", choices=cat_opts, selected=cat_opts)

    @reactive.effect
    @reactive.event(input.humu_gex_dot_splitby)  # Trigger when category changes
    def humu_update_dotplot_splitby():
        split_cat = input.humu_gex_dot_splitby()
        if split_cat != "---":
            cat_opts = hsh.categorial_opts_dict[split_cat]
            ui.update_selectize("humu_gex_dot_splits", choices=cat_opts, selected=cat_opts)
        else:
            ui.update_selectize("humu_gex_dot_splits", choices=[],)
    
    @render.plot
    @reactive.event(input.humu_gex_dot_run)
    def humu_plot_gex_dotplot():

        genes = list(input.humu_gex_dot_gene())
        if len(genes) > 0:
            groupby = input.humu_gex_dot_groupby()
            groups = input.humu_gex_dot_groups()
            splitby = input.humu_gex_dot_splitby()
            splits = input.humu_gex_dot_splits()
            swap = input.humu_gex_dot_swap()

            fig = hxp.plot_sc_dotplot(genes, groupby, groups, splitby, splits, swap)
            
            return fig
        else:
            return None
        
    
    @render_widget
    @reactive.event(input.humu_box_comp_run)
    def humu_gene_comparison_human():
        human_gene = input.humu_box_comp_human_genes()
        human_x_filter = input.humu_box_comp_x_cat_filter()
        human_transform = input.humu_box_comp_transformation()
        
        fig = hxp.humu_box_comparison_human(human_gene, human_x_filter, human_transform)

        return fig
    
    @render_widget
    @reactive.event(input.humu_box_comp_run)
    def humu_gene_comparison_mouse():
        mouse_gene = input.humu_box_comp_mu_genes()
        mouse_x_cat_filter = input.humu_box_comp_cat_mouse_subset()

        fig = hxp.humu_box_comparison_mouse(mouse_gene, mouse_x_cat_filter)

        return fig

    
    #@reactive.effect
    #@reactive.event(input.humu_box_comp_x_cat)  # Trigger when category changes
    #def update_box_viol_selectize():
    #    x_cat = input.humu_box_comp_x_cat()
    #    new_options = qsh.categorials_opts_dict[x_cat]
    #    ui.update_selectize("humu_box_comp_x_cat_filter", choices=new_options, selected=new_options)

        
    @reactive.effect
    @reactive.event(input.humu_box_comp_x_cat_mouse)
    def humu_update_box_viol_selectize():
        x_cat = input.humu_box_comp_x_cat_mouse()
        cat_opts = list(pd.read_feather("./quipi_humu_data/quipi_humu_adata_clean_full_PROC.feather", columns=[x_cat])[x_cat].unique())
        ui.update_selectize("humu_box_comp_cat_mouse_subset", choices=cat_opts, selected=cat_opts)
    



    



    @reactive.effect
    def populate_gene_selections():
        current_tab = input.humu_top_nav()
        if current_tab not in humu_tabs_visited:
            humu_tabs_visited.append(current_tab)
            if current_tab in humu_tabs_mapped_to_gene_inputs:
                for id, category in humu_tabs_mapped_to_gene_inputs[current_tab]:
                    if category == "Mouse":
                        ui.update_selectize(
                            id,
                            choices=hsh.genes,
                            selected=[],
                            server=True,
                        )
                    else:
                        ui.update_selectize(
                            id,
                            choices=hsh.quipi_genes,
                            selected=[],
                            server=True,
                        )

# Create the Shiny app
app_dir = Path(__file__).parent
app = App(app_ui, server,static_assets= app_dir / "www")



