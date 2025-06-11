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

gear_fill = ui.HTML(
    '<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" class="bi bi-gear-fill" viewBox="0 0 16 16"><path d="M9.405 1.05c-.413-1.4-2.397-1.4-2.81 0l-.1.34a1.464 1.464 0 0 1-2.105.872l-.31-.17c-1.283-.698-2.686.705-1.987 1.987l.169.311c.446.82.023 1.841-.872 2.105l-.34.1c-1.4.413-1.4 2.397 0 2.81l.34.1a1.464 1.464 0 0 1 .872 2.105l-.17.31c-.698 1.283.705 2.686 1.987 1.987l.311-.169a1.464 1.464 0 0 1 2.105.872l.1.34c.413 1.4 2.397 1.4 2.81 0l.1-.34a1.464 1.464 0 0 1 2.105-.872l.31.17c1.283.698 2.686-.705 1.987-1.987l-.169-.311a1.464 1.464 0 0 1 .872-2.105l.34-.1c1.4-.413 1.4-2.397 0-2.81l-.34-.1a1.464 1.464 0 0 1-.872-2.105l.17-.31c.698-1.283-.705-2.686-1.987-1.987l-.311.169a1.464 1.464 0 0 1-2.105-.872l-.1-.34zM8 10.93a2.929 2.929 0 1 1 0-5.86 2.929 2.929 0 0 1 0 5.858z"/></svg>'
)

humu_tabs_mapped_to_gene_inputs = {"Boxplots" : [["humu_gex_box_gene", "Mouse"]],
                                   "Dotplots" : [["humu_gex_dot_gene", "Mouse"]],
                                   "HuMu Expression Comparison" : [["humu_box_comp_human_genes", "Human"],["humu_box_comp_mu_genes", "Mouse"]],
}


#RUN_STYLE="background-color: #AFE1AF; color: black;"

panel_color = "#f0f0f0"


# Define the UI
app_ui = ui.page_fluid(

    ui.tags.style("""
        body { background-color: #8CCDEB; }  /* Background color */
        .nav-link { font-size: 20px;
                    color: white;
        }
        .navbar-nav .nav-link.active {
            color: white !important;  /* Keep text white when selected */
        }
                  
        .sidebar { /* This is a common class name for Shiny sidebars */
                  background-color: white !important; /* !important ensures it overrides other rules */
        }
    """),
    ui.head_content(
        ui.tags.link(
            rel="icon", type="image/png", sizes="64x64", href="favicon-32x32.png"
        ),
    ),
    ui.tags.div(
        ui.tags.div(
            ui.tags.img(src="humu.png", style="width: 90px; margin-left: 10px; margin-top: 18px; margin-bottom: 8px; border: 3px solid; border-radius: 8px;",class_="my-image"),  # Left-aligned image
            ui.tags.span("QuIPI" + "\n" + "HuMu", style="font-size: 35px; font-weight: bold;"),
            style="display: flex; align-items: flex-end; gap: 10px;" # Flexbox for horizontal alignment
        ),
        style="padding-bottom: 5px;"
    ),
    ui.page_navbar(

        ui.nav_panel("Home",
            ui.layout_column_wrap(
                ui.value_box(
                    "What is HuMu?",
                    "The Human-to-Mouse Cancer Translator Project (HuMu) aims to immuno-profile a series of \
                                    common and exceptional models of cancer in mice to benchmark them against the diversity of \
                                    TMEs in Human cancer (i.e immune archetypes described in Combes, Samad, et al. Cell 2022).",
                    theme=ui.value_box_theme(bg="#0B1D51" , fg="#FFFFFF"),
                    showcase_layout="top right",
                ),

                ui.value_box(
                    "HuMu Archetypes Translation Dataset",
                    ui.tags.img(src="quipi_humu_reference.png",
                                style="height: auto; border: 2px solid #ddd; border-radius: 10px; box-shadow: 5px 5px 15px rgba(0,0,0,0.3);"),
                                #theme="bg-gradient-teal-yellow",
                                theme=ui.value_box_theme(bg="#FFE3A9" , fg="#333333")
                ),
                ui.value_box(
                    "HuMu Dataset",
                    "The HuMu dataset contains high-dimensional cytometry data (CyTOF) of 15 murine models that \
                        we use to study high level tumor-immune composition, as well as single-cell sequencing data \
                        from 9 of these models that we use to dissect more granular gene expression profiles across \
                        populations and tumor models.",
                    theme=ui.value_box_theme(bg="#725CAD" , fg="#FFFFFF"),
                    showcase_layout="top right",
                ),
            )
        ),
        ui.nav_panel("HuMu Expression Comparison",
            ui.layout_sidebar(
                ui.sidebar(
                    ui.card(
                        ui.card_header("Human Options"),
                        ui.input_selectize("humu_box_comp_human_genes", "Choose Human Gene", [], multiple=False),
                        ui.input_selectize("humu_box_comp_x_cat_filter",
                                            "Subset Compartments",
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
                        ui.card_header("Mouse Options"),
                        ui.input_selectize("humu_box_comp_mu_genes", "Choose Murine Gene", []),
                        ui.input_selectize("humu_box_comp_cat_mouse_subset", "Subset Compartments",hsh.humu_compartments, selected=hsh.humu_compartments, multiple=True, remove_button=True,options={"plugins": ["clear_button"]}),
                        ui.input_switch("humu_box_comp_sample_aggr", "Average Counts by Sample")
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
            ),
        ),

        ui.nav_menu("Mouse Gene Expression",

            ui.nav_panel("Boxplots",
                ui.layout_sidebar(
                    ui.sidebar(
                        ui.input_selectize("humu_gex_box_gene", "Choose Gene to plot:", []),
                        ui.card(
                            ui.popover(
                                ui.span(
                                    gear_fill,
                                    style="position:absolute; top: 5px; right: 7px;",
                                ),
                                ui.input_selectize("humu_gex_box_cat_subset", "Subset Categories:", [], multiple=True, remove_button=True,options={"plugins": ["clear_button"]}),
                                placement="right",
                            ),
                            ui.input_selectize("humu_gex_box_x_cat", "Choose X-Axis Category:", hsh.categoricals_opts),
                        ),
                        ui.card(
                            ui.input_selectize("humu_gex_box_groupby", "Group by:", ["---"] + hsh.categoricals_opts, selected="---"),
                            ui.input_selectize("humu_gex_box_splitby", "Split by:", ["---"] + hsh.categoricals_opts, selected="---"),
                        ),
                        ui.card(ui.input_switch("humu_gex_box_sample_aggr", "Average Counts by Sample")),
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
                        ui.card(
                            ui.popover(
                                ui.span(
                                    gear_fill,
                                    style="position:absolute; top: 5px; right: 7px;",
                                ),
                                ui.input_selectize("humu_gex_dot_groups", "Subset Groupby Categories:", [], multiple=True),
                                placement="right",
                            ),
                            ui.input_selectize("humu_gex_dot_groupby", "Group By:", hsh.categoricals_opts, selected="---"),
                        ),
                        ui.card(
                            ui.popover(
                                ui.span(
                                    gear_fill,
                                    style="position:absolute; top: 5px; right: 7px;",
                                ),
                                ui.input_selectize('humu_gex_dot_splits', "Subset Splitby Categories:", [], multiple=True),
                                placement="right",
                            ),
                            ui.input_selectize("humu_gex_dot_splitby", "Split By:", ["---"] + hsh.categoricals_opts, selected="---"),
                        ),
                        ui.card(ui.input_switch("humu_gex_dot_swap", "Swap Axes")),
                        ui.input_action_button("humu_gex_dot_run", "RUN"),
                        bg=panel_color
                    ),
                ui.card(ui.output_plot("humu_plot_gex_dotplot"), full_screen=True),
                bg = panel_color
                ),       
            ),
        ),

        ui.nav_panel("Flow Cytometry Proportions",
            ui.layout_sidebar(
                ui.sidebar(
                    ui.card(
                        ui.input_selectize("humu_box_score_1", "Choose Parameter", hsh.flow_scores),
                        ui.input_selectize("humu_box_score_2", "Divide By Second Parameter (Optional)", ["---"] + hsh.flow_scores, selected="---"),
                    ),
                    ui.card(
                        ui.input_selectize("humu_box_x_cat", "Select X-Axis Category", hsh.flow_cats),
                        ui.input_switch("humu_box_switch", "Sort X-Category On Value"),
                            ui.popover(
                                ui.span(
                                    gear_fill,
                                    style="position:absolute; top: 5px; right: 7px;",
                                ),
                                ui.input_selectize("humu_box_x_cat_filter", "Subset X-Axis Categories.", [], multiple=True),
                                placement="right",
                            ),
                    ),
                    ui.card(ui.input_selectize("humu_box_group", "Group By:",  ["---"] + hsh.flow_cats, selected="---")),
                    ui.input_action_button("humu_box_run", "RUN"),
                    bg=panel_color
                ),
                ui.card(ui.card_body(output_widget("humu_expression_box_viol")),
                        ui.card_footer("Click button in the bottom right for fullscreen view."),
                        full_screen=True),
                bg=panel_color,
            )
        ),



        ##### EXAMPLE
        ui.nav_menu("Examples - Just Hit Run",
            ui.nav_panel("HuMu Expression Comparison", # EXAMPLE
                ui.layout_sidebar(
                    ui.sidebar(
                        ui.input_action_button("humu_box_comp_run_example", "RUN"),
                        ui.card(
                            ui.card_header("Human Options"),
                            ui.input_selectize("humu_box_comp_human_genes_example", "Choose Human Gene", ["CTLA4", "CCR5"], multiple=False, selected="CTLA4"),
                            ui.input_selectize("humu_box_comp_x_cat_filter_example",
                                                "Subset Compartments",
                                                hsh.quipi_compartments,
                                                selected=hsh.humu_compartments,
                                                multiple=True,
                                                remove_button=True,options={"plugins": ["clear_button"]}),
                            ui.input_selectize("humu_box_comp_transformation_example",
                                                "Select TPM Transformation:",
                                                ["TPM", "Log2(TPM)"],
                                                multiple= False,
                                                selected= "Log2(TPM)"),
                        ),
                        ui.card(
                            ui.card_header("Mouse Options"),
                            ui.input_selectize("humu_box_comp_mu_genes_example", "Choose Murine Gene", ["Ctla4", "Ccr5"], selected="Ctla4"),
                            ui.input_selectize("humu_box_comp_cat_mouse_subset_example", "Subset Categories:",hsh.humu_compartments, selected=hsh.humu_compartments, multiple=True, remove_button=True,options={"plugins": ["clear_button"]}),
                            ui.input_switch("humu_box_comp_sample_aggr_example", "Average Counts by Sample")
                        ),
                        bg=panel_color
                    ),
                    ui.card(
                        ui.card(ui.card_header("Human"),
                                output_widget("humu_gene_comparison_human_example")),
                        ui.card(ui.card_header("Mouse"),
                                output_widget("humu_gene_comparison_mouse_example")),
                        full_screen=True
                    ),
                bg=panel_color
                )
            ),
            ui.nav_panel("Mouse Boxplot", # EXAMPLE
                ui.layout_sidebar(
                    ui.sidebar(
                        ui.input_action_button("humu_gene_box_run_example", "RUN"),
                        ui.input_selectize("humu_gex_box_gene_example", "Choose Gene to plot:", ["Ctla4", "Ccr5"], selected="Ctla4"),
                        ui.card(
                            ui.input_selectize("humu_gex_box_x_cat_example", "Choose X-Axis Category:", hsh.categoricals_opts, selected = "Tumor Line"),
                            ui.popover(
                                    ui.span(
                                        gear_fill,
                                        style="position:absolute; top: 5px; right: 7px;",
                                    ),
                                    ui.input_selectize("humu_gex_box_cat_subset_example", "Subset Categories:", [], multiple=True, remove_button=True,options={"plugins": ["clear_button"]}),
                                    placement="right",
                                ),
                        ),
                        ui.card(
                            ui.input_selectize("humu_gex_box_groupby_example", "Group by:", ["---"] + hsh.categoricals_opts, selected="Coarse Annotation"),
                            ui.input_selectize("humu_gex_box_splitby_example", "Split by:", ["---"] + hsh.categoricals_opts, selected="---"),
                        ),
                        ui.card(ui.input_switch("humu_gex_box_sample_aggr_example", "Average Counts by Sample")),
                        bg=panel_color
                        ),
                        ui.card(output_widget("gex_box_example"), full_screen=True),
                        bg=panel_color
                )
            ),
            ui.nav_panel("Mouse Dotplot", # EXAMPLE
                ui.layout_sidebar(
                    ui.sidebar(
                        ui.input_selectize("humu_gex_dot_gene_example", "Choose Genes to plot:", 
                                           ["Itgam", "C1qc", "Cd3e", "Ctla4", "Pdcd1", "Lag3", "Havcr2", "Tox"], 
                                           selected=["Itgam", "C1qc", "Cd3e", "Ctla4", "Pdcd1", "Lag3", "Havcr2", "Tox"],multiple=True),
                        ui.card(
                            ui.popover(
                                ui.span(
                                    gear_fill,
                                    style="position:absolute; top: 5px; right: 7px;",
                                ),
                                ui.input_selectize("humu_gex_dot_groups_example", "Subset Groupby Categories:", [], multiple=True),
                                placement="right",
                            ),
                            ui.input_selectize("humu_gex_dot_groupby_example", "Group By:", hsh.categoricals_opts, selected="Compartment"),
                        ),
                        ui.card(
                            ui.popover(
                                ui.span(
                                    gear_fill,
                                    style="position:absolute; top: 5px; right: 7px;",
                                ),
                                ui.input_selectize('humu_gex_dot_splits_example', "Subset Splitby Categories:", [], multiple=True),
                                placement="right",
                            ),
                            ui.input_selectize("humu_gex_dot_splitby_example", "Split By:", ["---"] + hsh.categoricals_opts, selected="Tumor Line"),
                        ),
                        ui.card(ui.input_switch("humu_gex_dot_swap_example", "Swap Axes", True)),
                        ui.input_action_button("humu_gex_dot_run_example", "RUN"),
                        bg=panel_color
                    ),
                ui.card(ui.output_plot("humu_plot_gex_dotplot_example"), full_screen=True),
                bg = panel_color
                )
            ),
            ui.nav_panel("Flow Cytometry Proportions", # EXAMPLE
                ui.layout_sidebar(
                    ui.sidebar(
                        ui.input_action_button("humu_box_run_example", "RUN"),
                        ui.card(
                            ui.input_selectize("humu_box_score_1_example", "Choose Parameter", ["Tcells_Live"], selected = "Tcells_Live"),
                            ui.input_selectize("humu_box_score_2_example", "Divide By Second Parameter (Optional)", ["CD45_Live"], selected="CD45_Live"),
                        ),
                        ui.card(
                            ui.input_selectize("humu_box_x_cat_example", "Select X-Axis Category", hsh.flow_cats, selected="Human Archetype / Mouse Tumor Line"),
                            ui.popover(
                                    ui.span(
                                        gear_fill,
                                        style="position:absolute; top: 5px; right: 7px;",
                                    ),
                                    ui.input_selectize("humu_box_x_cat_filter_example", "Subset X-Axis Categories.", [], multiple=True),
                                    placement="right",
                                ),
                            ui.input_switch("humu_box_switch_example", "Sort X-Category On Value", True),
                        ),
                        ui.card(ui.input_selectize("humu_box_group_example", "Group By:",  ["---"] + hsh.flow_cats, selected="Species")),
                        bg=panel_color
                    ),
                ui.card(ui.card_body(output_widget("humu_expression_box_viol_example")),
                        ui.card_footer("Click button in the bottom right for fullscreen view."),
                        full_screen=True),
                bg=panel_color,
            )
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
    
    ##### GEX Box Plots
    @render_widget
    @reactive.event(input.humu_gene_box_run)
    def gex_box():
        gene = input.humu_gex_box_gene()
        x_cat = input.humu_gex_box_x_cat()
        x_sub = input.humu_gex_box_cat_subset()
        group = input.humu_gex_box_groupby()
        split = input.humu_gex_box_splitby()
        sample_aggr = input.humu_gex_box_sample_aggr()

        fig = hxp.plot_sc_box(gene,x_cat,x_sub,group,split, sample_aggr)
        return fig

    @reactive.effect
    @reactive.event(input.humu_gex_box_x_cat)  # Trigger when category changes
    def humu_update_box_viol_selectize():
        x_cat = input.humu_gex_box_x_cat()
        cat_opts = hsh.categorial_opts_dict[x_cat]
        #cat_opts = list(pd.read_feather("./quipi_humu_data/quipi_humu_adata_clean_full_PROC.feather", columns=[x_cat])[x_cat].unique())
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
        

    ##### HUMU COMPARISON BOXPLOTS   
    
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
        mouse_sample_aggr = input.humu_box_comp_sample_aggr()

        fig = hxp.humu_box_comparison_mouse(mouse_gene, mouse_x_cat_filter, mouse_sample_aggr)

        return fig

        
    @reactive.effect
    @reactive.event(input.humu_box_comp_x_cat_mouse)
    def humu_update_box_viol_selectize():
        x_cat = input.humu_box_comp_x_cat_mouse()
        cat_opts = list(pd.read_feather("./quipi_humu_data/quipi_humu_adata_clean_full_PROC.feather", columns=[x_cat])[x_cat].unique())
        ui.update_selectize("humu_box_comp_cat_mouse_subset", choices=cat_opts, selected=cat_opts)



    ##### EXAMPLES


    ## COMPARISON
    @render_widget
    @reactive.event(input.humu_box_comp_run_example)
    def humu_gene_comparison_human_example():
        human_gene = input.humu_box_comp_human_genes_example()
        human_x_filter = input.humu_box_comp_x_cat_filter_example()
        human_transform = input.humu_box_comp_transformation_example()


        fig = hxp.humu_box_comparison_human(human_gene, human_x_filter, human_transform)

        return fig
    
    @render_widget
    @reactive.event(input.humu_box_comp_run_example)
    def humu_gene_comparison_mouse_example():
        mouse_gene = input.humu_box_comp_mu_genes_example()
        mouse_x_cat_filter = input.humu_box_comp_cat_mouse_subset_example()
        mouse_sample_aggr = input.humu_box_comp_sample_aggr_example()

        fig = hxp.humu_box_comparison_mouse(mouse_gene, mouse_x_cat_filter, mouse_sample_aggr)

        return fig
    
    ## BOXPLOTS
    @render_widget
    @reactive.event(input.humu_gene_box_run_example)
    def gex_box_example():
        gene = input.humu_gex_box_gene_example()
        x_cat = input.humu_gex_box_x_cat_example()
        x_sub = input.humu_gex_box_cat_subset_example()
        group = input.humu_gex_box_groupby_example()
        split = input.humu_gex_box_splitby_example()
        sample_aggr = input.humu_gex_box_sample_aggr_example()

        fig = hxp.plot_sc_box(gene,x_cat,x_sub,group,split, sample_aggr)
        return fig

    @reactive.effect
    @reactive.event(input.humu_gex_box_x_cat_example)  # Trigger when category changes
    def humu_update_box_viol_selectize():
        x_cat = input.humu_gex_box_x_cat_example()
        cat_opts = hsh.categorial_opts_dict[x_cat]
        
        #list(pd.read_feather("./quipi_humu_data/quipi_humu_adata_clean_full_PROC.feather", columns=[x_cat])[x_cat].unique())
        ui.update_selectize("humu_gex_box_cat_subset_example", choices=cat_opts, selected=cat_opts)

    ## FLOWPLOTS

    @reactive.effect
    @reactive.event(input.humu_box_x_cat_example)  # Trigger when category changes
    def update_box_viol_selectize():
        x_cat = input.humu_box_x_cat_example()
        flow_df = pd.read_feather("./quipi_humu_data/quipi_humu_flow_table.feather", columns=[x_cat])
        cats = list(flow_df[x_cat].unique())
        ui.update_selectize("humu_box_x_cat_filter_example", choices=cats, selected=cats)

    @render_widget
    @reactive.event(input.humu_box_run_example)
    def humu_expression_box_viol_example():
        score1 = input.humu_box_score_1_example()
        score2 = input.humu_box_score_2_example()
        x_cat = input.humu_box_x_cat_example()
        x_cat_filter = input.humu_box_x_cat_filter_example()
        group = input.humu_box_group_example()
        sort = input.humu_box_switch_example()
        fig = hfb.box_humu_flow(score1, score2, x_cat, x_cat_filter, group, sort)
        return fig
    
    ## DOTPLOTS
    @reactive.effect
    @reactive.event(input.humu_gex_dot_groupby_example)
    def humu_update_dotplot_groupby_example():
        x_cat = input.humu_gex_dot_groupby_example()
        cat_opts = hsh.categorial_opts_dict[x_cat]
        ui.update_selectize("humu_gex_dot_groups_example", choices=cat_opts, selected=cat_opts)

    @reactive.effect
    @reactive.event(input.humu_gex_dot_splitby_example)  # Trigger when category changes
    def humu_update_dotplot_splitby_example():
        split_cat = input.humu_gex_dot_splitby_example()
        if split_cat != "---":
            cat_opts = hsh.categorial_opts_dict[split_cat]
            ui.update_selectize("humu_gex_dot_splits_example", choices=cat_opts, selected=cat_opts)
        else:
            ui.update_selectize("humu_gex_dot_splits_example", choices=[],)
    
    @render.plot
    @reactive.event(input.humu_gex_dot_run_example)
    def humu_plot_gex_dotplot_example():

        genes = list(input.humu_gex_dot_gene_example())
        if len(genes) > 0:
            groupby = input.humu_gex_dot_groupby_example()
            groups = input.humu_gex_dot_groups_example()
            splitby = input.humu_gex_dot_splitby_example()
            splits = input.humu_gex_dot_splits_example()
            swap = input.humu_gex_dot_swap_example()

            fig = hxp.plot_sc_dotplot(genes, groupby, groups, splitby, splits, swap)
            
            return fig
        else:
            return None



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



