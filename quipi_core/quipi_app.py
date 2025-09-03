from shiny import App, render, ui, reactive
from shinyswatch import theme
import plotly.express as px

px.defaults.template = "simple_white"
from shinywidgets import output_widget, render_widget 

import quipi_shared as qsh
#import humu_shared as hsh
import pandas as pd
from io import StringIO

import ranked_patient_dge as rpd
import gene_factor as gf
import box_viol_expression_plot as bv
import pancan_plots as pp
import correlation_analysis as corr
from pathlib import Path

from faicons import icon_svg

#import humu_flow_boxplot as hfb
#import humu_gex_plots as hxp

gear_fill = ui.HTML(
    '<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" class="bi bi-gear-fill" viewBox="0 0 16 16"><path d="M9.405 1.05c-.413-1.4-2.397-1.4-2.81 0l-.1.34a1.464 1.464 0 0 1-2.105.872l-.31-.17c-1.283-.698-2.686.705-1.987 1.987l.169.311c.446.82.023 1.841-.872 2.105l-.34.1c-1.4.413-1.4 2.397 0 2.81l.34.1a1.464 1.464 0 0 1 .872 2.105l-.17.31c-.698 1.283.705 2.686 1.987 1.987l.311-.169a1.464 1.464 0 0 1 2.105.872l.1.34c.413 1.4 2.397 1.4 2.81 0l.1-.34a1.464 1.464 0 0 1 2.105-.872l.31.17c1.283.698 2.686-.705 1.987-1.987l-.169-.311a1.464 1.464 0 0 1 .872-2.105l.34-.1c1.4-.413 1.4-2.397 0-2.81l-.34-.1a1.464 1.464 0 0 1-.872-2.105l.17-.31c.698-1.283-.705-2.686-1.987-1.987l-.311.169a1.464 1.464 0 0 1-2.105-.872l-.1-.34zM8 10.93a2.929 2.929 0 1 1 0-5.86 2.929 2.929 0 0 1 0 5.858z"/></svg>'
)

quipi_tabs_mapped_to_gene_inputs = {"Boxplot" : ["box_viol_gene_input"],
                              "Gene Expresssion Bar Plots" : ["gene_expr_bar_genes"],
                              "Correlation Matrix" : ["corr_gene_input"],
                              "Compartment Correlation Matrix" : ["comp_corr_mat_genes"],
                              "One-Vs-All Correlation Table" : ["corr_cat_gene_input"],
                              "PanCan UMAP Gene Expression" : ["pancan_gene_input"],
                              "PanCan Gene-Signature Overlay" : ["gene_factor_genes"],
                              "Feature Score Ranked DGE" : ["dge_highlight_genes"],
                              "Gene-Signature Score Ranked DGE" : ["fs_dge_genes","fs_dge_highlight_genes"],
                              "Query Data" : ["gene_expr_query_genes"],
                              "Cross-Compartment Correlation": ["cross_comp_corr_cat_gene_input"],
                              "Ridgeline" : ["ridgeline_gene_input"]

}

panel_color = "#f0f0f0"
input_min_height = 200
#sidebar_color = "#8DC8AF"
#panel_color = "#BDD0D5"


# Define the UI
app_ui = ui.page_navbar(
        
    ## HOME PAGE
    ui.nav_panel( # QUIPI HOME PAGE
        ui.tags.img(src="quipi.png", style="height: 55px; display: inline-block; vertical-align: middle;"),
        ui.card(

            # QuIPI Description
            ui.card(
                    ui.tags.div(
                    ui.tags.img(src="quipi.png", style="height: 100px; margin-left: 10px; padding-top: 0px"),  # Left-aligned image
                ),
                ui.h2("QuIPI - Querying Immunoprofiler"),
                ui.h4(
                    """
                    QuIPI is a user interface designed to make it easy to visualize the UCSF Immunoprofiler dataset. 
                    QuIPI incorporates Bulk RNA-Seq and Flow Cytometry data from over 400 cancer samples across multiple
                    indications. QuIPI also incorporates archetype calls, building off of the work of Combes, Samad et al. Cell 2022 
                    in order to identify common immune states across cancer types.

                    QuIPI provides quick and easy to use plotting functionality, allowing anyone to discover trends and insights from
                    Immunoprofiler data. QuIPI modules are designed with the hope that they can be generalizable to different questions
                    while still allowing for sufficient cusomization.
                    """
                ),
                ui.HTML(f'<a href="{"https://pubmed.ncbi.nlm.nih.gov/34963056/"}" target="_blank">Discovering dominant tumor immune archetypes in a pan-cancer census. Combes AJ, Samad B, Tsui J, et al. Cell. 2022</a>'),
                style="background-color: #BDD0D5"
            ),
            # Immunoprofiler Description
            ui.card(
                ui.tags.div(
                    ui.tags.img(src="immunoprofiler_image.png", style="height: 100px; margin-left: 10px; padding-top: 10px"),  # Left-aligned image
                    ui.tags.span("", style="font-size: 50px; font-weight: bold;"),
                    style="display: flex; align-items: flex-end; gap: 10px;" # Flexbox for horizontal alignment
                ),
                ui.h2("Immunoprofiler"),
                ui.h4(
                    """
                    The Immunoprofiler (IPI) effort coordinates the processing of valuable samples from cancer patients
                    in order to perform various tests to learn more about immune composition, immune cell gene expression, and immune interaction biology.
                    By viewing cancers as a discrete forms of immunopathology, better understanding the immune response is a critical piece for
                    how to treat cancers and provide new targets for the new forms of immunotherapies.
                    """
                ),
                ui.HTML(f'<a href="{"https://immunoprofiler.org/cancer"}" target="_blank">UCSF Immunoprofiler</a>'),
                style="background-color: #8DC8AF"
            ),
        )
    ),

    # ABOUT Panel
    ui.nav_panel("About The Dataset",
                 
        ui.accordion(

            # Plots describing cancer indications
            ui.accordion_panel("Cancer Indications",
                ui.layout_column_wrap(
                    ui.card(ui.card_header("Cancer Indication Breakdown"),
                            output_widget("plot_indication_breakdown"), 
                            full_screen=False,
                            ),
                    ui.card(ui.card_header("Cancer Indication Nomenclature"),
                        output_widget("cancer_glossary"),
                        full_screen=False),
                        width=.5, min_height=650
                ),
            ),
            # Plots describing pancan archetypes
            ui.accordion_panel("PanCan Archetypes",
                ui.layout_column_wrap(
                    ui.card(ui.card_header("PanCan Archetype UMAP"),
                            ui.card_body(output_widget("pancan_archetypes_home")),
                            ui.card_footer(
                                ui.HTML(f'<a href="{"https://pubmed.ncbi.nlm.nih.gov/34963056/"}" target="_blank">Discovering dominant tumor immune archetypes in a pan-cancer census. Combes AJ, Samad B, Tsui J, et al. Cell. 2022</a>')
                            ),
                            full_screen=False,     
                    ),
                    ui.card(ui.card_header("PanCan Archetype Breakdown"),
                            ui.card_body(output_widget("pancan_archetype_breakdown")),
                    ),
                    min_height=700
                ),
            ),
            id="data_accordion",
            open=False
        ),
    ),

    ##### GENE EXPRESSION PLOTS

    # BOXPLOT NAV
    ui.nav_panel("Boxplot",
        ui.layout_sidebar(
            ui.sidebar(
                ui.h3("Boxplot"),
                ui.accordion(
                    ui.accordion_panel("Select Genes",
                        ui.input_selectize("box_viol_gene_input",
                                            "Select Genes:",
                                            [],
                                            multiple=True),
                        ui.output_ui("box_viol_multiple_genes"),
                    ),
                    ui.accordion_panel("X-Axis Category",
                        ui.popover(
                            ui.span(gear_fill, style="float: right;"),
                            ui.input_selectize("box_viol_x_cat_filter",
                                    "Subset X-Axis Categories",
                                    [],
                                    multiple=True,
                                    remove_button=True,options={"plugins": ["clear_button"]}),
                            placement="right",
                        ),
                        ui.input_selectize("box_viol_x_category",
                                        "",
                                        list(qsh.categoricals_dict.keys()),
                                        selected = "Compartment"),
                    ),

                    ui.accordion_panel("Groupby:",
                        ui.input_selectize("box_viol_groupby",
                            "",
                            ["---"] + list(qsh.categoricals_dict.keys()),
                            selected="---")),

                    ui.accordion_panel("TPM Transformation",
                        ui.input_selectize("box_viol_transformation",
                            "",
                            ["TPM", "Log2(TPM)"],
                            multiple= False,
                            selected= "Log2(TPM)")),
                ),
                    ui.input_action_button("box_viol_run", "RUN", icon=icon_svg("arrow-right")),
                    bg=panel_color
            ),
            ui.card((output_widget("expression_box_viol")),full_screen=True),
            bg=panel_color
        ),
    ),


    # # RIDGELINE PLOT NAV
    # ui.nav_panel("Ridgeline",
    #     ui.layout_sidebar(
    #         ui.sidebar(
    #             ui.h4("Ridgeline Plot"),
    #             ui.card(
    #                 ui.input_selectize("ridgeline_gene_input",
    #                                     "Select Genes:",
    #                                     [],
    #                                     multiple=True)),
    #                 ui.card(
    #                     ui.popover(
    #                         ui.span(gear_fill, style="position:absolute; top: 5px; right: 7px;"),
    #                         ui.input_selectize("ridgeline_x_cat_filter",
    #                                 "**Subset X-Axis Categories.**",
    #                                 [],
    #                                 multiple=True,
    #                                 remove_button=True,options={"plugins": ["clear_button"]}),
    #                         placement="right",
    #                     ),
    #                     ui.input_selectize("ridgeline_x_category",
    #                                     "Select X-Axis Category:",
    #                                     list(qsh.categoricals_dict.keys()),
    #                                     selected = "Compartment"),
    #                 ),
    #                 ui.card(ui.input_selectize("ridgeline_groupby",
    #                                     "Group by:",
    #                                     ["---"] + list(qsh.categoricals_dict.keys()),
    #                                     selected="---")),
    #                 ui.card(ui.input_selectize("ridgeline_transformation",
    #                                     "TPM Transformation: ",
    #                                     ["TPM", "Log2(TPM)"],
    #                                     multiple= False,
    #                                     selected= "Log2(TPM)")),
    #                 ui.input_action_button("ridgeline_run", "RUN",icon=icon_svg("arrow-right")),
    #                 bg=panel_color
    #         ),
    #     ui.card(ui.output_ui("ridgeline"), full_screen=True),
    #     bg=panel_color,
    #     ),       
    # ),

    # HEATMAP NAV
    ui.nav_panel("Heatmap",  
        ui.layout_sidebar(
            ui.sidebar(
                ui.h4("Heatmap"),

                ui.accordion(
                    ui.accordion_panel("Input Genes",
                        ui.input_text_area(
                        "gex_heatmap_text_genes",
                        "",
                        rows=8
                        ),
                    ),
                    ui.accordion_panel("Choose X-Axis Category",
                        ui.popover(
                            ui.span(gear_fill,style="float: right;"),
                            ui.input_selectize("gex_heatmap_groupby_subset", "Subset Categories:", [], multiple=True,remove_button=True,options={"plugins": ["clear_button"]}),
                            placement="right",
                        ),
                        ui.input_selectize("gex_heatmap_groupby", "",list(qsh.categorical_choices.keys())),
                    ),
                    ui.accordion_panel("TPM Transformation",   
                        ui.input_selectize("gex_heatmap_transform", "", ["Log2(TPM)", "TPM"], selected = "TPM")
                    ),
                ),
                ui.input_action_button("gex_heatmap_run", "RUN", icon=icon_svg("arrow-right")),
                bg=panel_color
            ),
            ui.card(ui.output_ui("gex_heatmap"), full_screen=True),
            bg=panel_color,
        ),
    ),

    # DATA QUERY NAV
    ui.nav_panel("Data Query",
        ui.layout_sidebar(
            ui.sidebar(
                ui.h4("Data Query"),
                ui.accordion(
                    ui.accordion_panel("Input Genes",
                        ui.input_text_area(
                            "gene_expr_query_genes",
                            "",
                            rows=4
                        ),
                    ),
                    ui.accordion_panel("Select Indications",
                        ui.input_selectize("query_indication",
                                            "",
                                            choices=qsh.indications,
                                            selected=qsh.indications,
                                            multiple = True,
                                            remove_button=True,options={"plugins": ["clear_button"]}
                        ),
                    ),
                    ui.accordion_panel("Select Compartments",
                        ui.input_selectize("query_compartment",
                                            "",
                                            choices=qsh.compartments,
                                            selected=qsh.compartments,
                                            multiple=True,
                                            remove_button=True,options={"plugins": ["clear_button"]}
                        ),
                    ),
                    ui.accordion_panel("Select Archetypes",
                        ui.input_selectize("query_archetype",
                                            "",
                                            choices=qsh.archetypes,
                                            multiple=True,
                                            selected=qsh.archetypes,
                                            remove_button=True,options={"plugins": ["clear_button"]}
                        ),
                    ),
                    ui.accordion_panel("TPM Transformation",
                        ui.input_selectize("query_transform",
                                            "",
                                            choices=["TPM", "Log2(TPM)"],
                                            selected="Log2(TPM)"
                        ),
                    ),
                ),
            ui.input_action_button("gene_expr_query_run", "RUN", icon=icon_svg("arrow-right")),
            ui.download_button("download_query_table", "Download CSV", icon=icon_svg("download")),
            bg=panel_color
            ),
            ui.card(ui.output_data_frame("gene_expr_query"),full_screen=True),
            bg=panel_color
        ),
    ),

    # CORRELATION MENU
    ui.nav_menu("Correlation", 
        
        # BASIC CORRELATION MATRIX
        ui.nav_panel("Correlation Matrix",
            ui.h4("Explore the correlation of genes across the entire IPI dataset."),
            ui.layout_sidebar(
                ui.sidebar(
                    ui.accordion(
                        ui.accordion_panel("Select Genes",
                            ui.input_selectize("corr_gene_input",
                                                "",
                                                [],
                                                multiple=True,
                                                options = {'server':True}
                            ),
                        ),
                        ui.accordion_panel("Select Indications",
                            ui.input_selectize("corr_indication",
                                            "",
                                            choices=qsh.indications,
                                            selected=qsh.indications,
                                            multiple = True,
                                            remove_button=True,options={"plugins": ["clear_button"]}
                            ),
                        ),
                        ui.accordion_panel("Select Tissue",
                            ui.input_selectize("corr_tissue",
                                                "",
                                                choices=["Tumor", "Normal"],
                                                selected = ["Tumor", "Normal"],
                                                multiple = True,
                                                remove_button=True,options={"plugins": ["clear_button"]}
                            ),
                        ),
                        ui.accordion_panel("Select Compartments",
                            ui.input_selectize("corr_compartment",
                                            "",
                                            choices=qsh.compartments,
                                            selected=qsh.compartments,
                                            multiple=True,
                                            remove_button=True,options={"plugins": ["clear_button"]}
                            ),
                        ),
                        ui.accordion_panel("Select Archetypes",
                            ui.input_selectize("corr_archetype",
                                                "",
                                                choices=qsh.archetypes,
                                                multiple=True,
                                                selected=qsh.archetypes,
                                                remove_button=True,options={"plugins": ["clear_button"]}
                            )
                        ),
                        ui.accordion_panel("TPM Transformation",
                            ui.input_selectize("corr_transform",
                                                "",
                                                choices=["TPM", "Log2(TPM)"],
                                                selected="Log2(TPM)"
                            ),
                        ),

                        ui.accordion_panel("Correlation Method",
                            ui.input_selectize("corr_method_input",
                                                "",
                                                choices=["Pearson", "Spearman"],
                                                selected="Spearman"
                            ),
                        )
                    ),
                    ui.input_action_button("corr_run", "RUN",icon=icon_svg("arrow-right")),
                    bg=panel_color
                ),
                ui.card(ui.card_body(output_widget("gene_correlation_heatmap")),full_screen=True),
                bg=panel_color
            )
        ),

        # # CROSS-COMPARTMENT CORRELATION MATRIX
        # ui.nav_panel("Compartment Correlation Matrix",
        #     ui.h4("Explore the correlation of chosen genes across two compartments."),
        #     ui.layout_sidebar(
        #         ui.sidebar(
        #             ui.input_selectize("comp_corr_mat_genes",
        #                             "Select Genes:",
        #                             [],
        #                             multiple=True,
        #                             options = {"server":True}),
        #             ui.input_selectize("comp_corr_mat_compartment1",
        #                                 "Select Compartment A:",
        #                                 choices=qsh.compartments,
        #                                 ),
        #             ui.input_selectize("comp_corr_mat_compartment2",
        #                             "Select Compartment B:",
        #                             choices=qsh.compartments),
        #             ui.input_selectize("comp_corr_mat_method",
        #                             "Select Correlation Method:",
        #                             choices = ["Pearson", "Spearman"],
        #                             selected="Spearman"),
        #             ui.input_selectize("comp_corr_mat_transform",
        #                                 "Select TPM Transformation:",
        #                                 choices=["TPM", "Log2(TPM)"],
        #                                 selected="Log2(TPM)"),
        #             ui.input_selectize("comp_corr_mat_indication",
        #                             "Select Indications:",
        #                             choices=qsh.indications,
        #                             selected=qsh.indications,
        #                             multiple = True),
        #             ui.input_selectize("comp_corr_mat_tissue",
        #                                 "Select Tissue:",
        #                                 choices=["Tumor", "Normal"],
        #                                 selected = ["Tumor", "Normal"],
        #                                 multiple = True),
        #             ui.input_selectize("comp_corr_mat_archetype",
        #                                 "Select Archetypes:",
        #                                 choices=qsh.archetypes,
        #                                 multiple=True,
        #                                 selected=qsh.archetypes),
        #             ui.input_action_button("comp_corr_mat_run", "RUN",icon=icon_svg("arrow-right")),
        #             bg=panel_color
        #         ),
        #         ui.card(ui.card_body(output_widget("comp_corr_mat")),full_screen=True),
        #         bg=panel_color
        #     )
        # ),
        
        # # CORRELATION OF 1 GENE WITH ALL OTHER GENES IN A PARTICULAR CATEGORY
        # ui.nav_panel("One-Vs-All Correlation Table",
        #     ui.h3("One-Vs-All correlation analysis within a chosen category."),
        #     ui.layout_sidebar(
        #             ui.sidebar(
        #                 ui.input_selectize("corr_cat_gene_input",
        #                                     "Select Gene:",
        #                                     [],
        #                                     multiple=False,
        #                                     options = {'server':True}),
        #                 ui.input_selectize("corr_cat_category_input",
        #                                 "Select Category:",
        #                                 choices=list(qsh.categorical_choices.keys())),
        #                 ui.input_selectize("corr_cat_category_opts",
        #                                 "Select Subcategories:",
        #                                 choices = [],
        #                                 multiple=True),
        #                 ui.input_slider("corr_cat_slider", "Select Correlation Coefficient Range", min=-1,max=1, value = (0,1), step=.1),
        #                 ui.input_action_button("corr_cat_run", "RUN",icon=icon_svg("arrow-right")),
        #                 bg=panel_color           
        #             ),
        #         ui.card(ui.output_data_frame("corr_cat_gene_correlations"), full_screen=True),
        #         bg=panel_color
        #     )
        # ),
        
        # CORRELATION OF GENES IN COMPARTMENT A WITH ALL GENES IN COMPARTMENT B
        ui.nav_panel("Cross-Compartment Correlation",
            ui.layout_sidebar(
                ui.sidebar(
                    ui.h5("Cross-Compartment Correlation"),
                    ui.accordion(
                        ui.accordion_panel("Gene",
                            ui.input_selectize("cross_comp_corr_cat_gene_input",
                                                "",
                                                [],
                                                multiple=True,
                                                options = {'server':True}
                            ),
                        ),
                        ui.accordion_panel("Compartment A",
                            ui.input_selectize("cross_comp_corr_comp1_input",
                                            "",
                                            choices=list(qsh.compartments)
                            ),
                        ),
                        ui.accordion_panel("Compartment B",
                            ui.input_selectize("cross_comp_corr_comp2_input",
                                            "",
                                            choices = list(qsh.compartments),
                                            multiple=True),
                        ),
                        ui.accordion_panel("Correlation Coefficient Range",
                            ui.input_slider("cross_comp_corr_cat_slider", "", 
                                            min=-1,max=1, value = (.5,1), step=.1
                            )
                        ),
                        ui.accordion_panel("Correlation Methods",
                            ui.input_selectize("cross_comp_corr_mat_method",
                                        "",
                                        choices = ["Pearson", "Spearman"],
                                        selected="Spearman"
                            )
                        ),
                        ui.accordion_panel("TPM Transformation",
                            ui.input_selectize("cross_comp_corr_mat_transform",
                                            "",
                                            choices=["TPM", "Log2(TPM)"],
                                            selected="Log2(TPM)"
                            ),
                        )
                    ),
                    ui.input_action_button("cross_comp_corr_cat_run", "RUN",icon=icon_svg("arrow-right")),
                    bg=panel_color
                ),
                ui.card(ui.output_data_frame("cross_comp_corr_cat_gene_correlations")),
                bg=panel_color
            )
        )
    ), # END CORRELATION MENU

    # PANCAN ARCHETYPES MENU
    ui.nav_menu("PanCan Archetypes",
        # Tab where user can select multiple genes and view their expression overlayed on the individual
        ui.nav_panel("PanCan UMAP Gene Expression",
            ui.layout_sidebar(
                ui.sidebar(
                    ui.h5("PanCan Gene Expression Overlay"),

                    ui.accordion(

                        ui.accordion_panel("Genes",
                            ui.input_selectize("pancan_gene_input",
                                                "",
                                                [],
                                                multiple = True
                            ),
                        ),
                        ui.accordion_panel("Compartment",
                            ui.input_selectize("pancan_compartment_input",
                                            "",
                                            qsh.compartments
                            ),
                        ),
                        ui.accordion_panel("TPM Transformation",
                            ui.input_selectize("pancan_umap_transformation",
                                            "",
                                            ["TPM", "Log2(TPM)"],
                                            multiple= False,
                                            selected= "Log2(TPM)"
                            ),
                        )
                    ),
                    ui.input_action_button("pancan_umap_run", "RUN",icon=icon_svg("arrow-right")),
                    bg=panel_color,
                ),
                ui.card(output_widget("pancan_subplots")),
                bg=panel_color
            ),
        ),

        # Given user-defined gene list, plots the gene factor score on top of the PanCan UMAP
        ui.nav_panel("PanCan Gene-Signature Overlay",
            ui.layout_sidebar(
                ui.sidebar(
                    ui.h5("PanCan Gene-Signature Score Overlay"),

                    ui.accordion(
                        ui.accordion_panel("Genes",
                            ui.input_selectize("gene_factor_genes",
                                                "",
                                                [],
                                                multiple=True
                            )
                        ),
                        ui.accordion_panel("Compartment",
                            ui.input_selectize("gene_factor_compartment",
                                            "",
                                            qsh.compartments
                            ),
                        )
                    ),
                    ui.input_action_button("gene_factor_run", "RUN",icon=icon_svg("arrow-right")),
                    bg=panel_color
                ),
                ui.layout_column_wrap(
                    ui.card(ui.card_body(output_widget("gene_factor_analysis"))),
                    ui.card(ui.card_body(output_widget("pancan_archetypes_gfs")))
                ),
                bg=panel_color
            ),
            ui.h3("Gene-Signature Score Calculation"),
            ui.p("""
                Data is subset to include archetyped patients only within a given compartment.
                Gene-signature scores are calculated by first calculating a per-gene Z-Score across patients
                and subsequently averaging the Z-Scores within each patient resulting in each patient being
                assigned a gene-signature score depending on their expression of the input gene set.

                """)
        )
    ), # END PANCAN MENU

    # BEGIN DGE MENU
    ui.nav_menu("Differential Gene Expression",

        # FEATURE SCORE RANKED DGE
        ui.nav_panel("Feature Score Ranked DGE",
            ui.layout_sidebar(
                ui.sidebar(
                    ui.h4("Feature Score Ranked DGE"),
                    ui.accordion(
                        ui.accordion_panel("Feature Score for ranking",
                            ui.input_selectize("flow_score_to_rank",
                                            "",
                                            list(qsh.feature_scores.keys())
                            ),
                        ),
                        ui.accordion_panel("Genes to highlight",
                            ui.input_selectize("dge_highlight_genes",
                                            "",
                                            [],
                                            multiple=True
                            ),
                        ),
                        ui.accordion_panel("Quartile",
                            ui.input_slider("dge_slider",
                                            "Selects the top and bottom Quartile% of patients",
                                            value = .2,
                                            min = 0.05,
                                            max = .5,
                                            step = .05
                            ),
                        ),
                        ui.accordion_panel("Compartment",
                            ui.input_selectize("dge_compartment",
                                                "DEGs calculated within this compartment",
                                                qsh.compartments
                            )
                        ),
                        ui.accordion_panel("Thresholds",
                            ui.input_numeric("dge_fc_thresh",
                                            "Fold Change Magnitude",
                                            value=2
                            ),
                            ui.input_numeric("dge_p_thresh",
                                            "-Log10(P-Value)",
                                            value = 1
                            ),
                        )
                    ),
                    ui.input_action_button("dge_run", "RUN",icon=icon_svg("arrow-right")),
                    bg=panel_color
                ),
                ui.card(
                    ui.card(output_widget("compartment_featurescore_dge")),
                    ui.layout_column_wrap(
                        ui.card(
                            ui.card_header("Negative DEGs"),
                            ui.card_body(ui.output_data_frame("compartment_featurescore_dge_bot")),
                        ),
                        ui.card(
                            ui.card_header("Positive DEGs"),
                            ui.card_body(ui.output_data_frame("compartment_featurescore_dge_top"))
                        ),
                    ),
                    ui.download_button("download_featurescore_dges", "Download DEGs", icon=icon_svg("download")),
                ),
                ui.card(
                    ui.card_header("Expression Levels for genes of interest"),
                    ui.card_body(output_widget("compartment_featurescore_dge_boxplot"))
                ),
                bg=panel_color
            ),

            # METHODS SECTION
            ui.h4("Feature Scoring and DGE Group Formation"),
            ui.p("""
                Data is subset to include only patients assigned feature scores. Effectively, this means that
                patients without an archetype designation are not considered. Patients are ranked based
                on their feature score in the specified compartment and are split based on the specified
                quartile. For the chose quartile, the resulting groups used for DGE are top quartile percent
                of patients with the highest score and the bottom quartile percent of patients with the lowest score.
            """),
            ui.h4("Differential Gene Expression"),
            ui.p("""The Wilcoxon Rank-sum test is used to calculate differentially expressed genes.
                    P-values are corrected using Benjamini/Hochberg multiple testing correction.
                    """),
        ),

        # DGE CALCULATION BASED ON TOP VS. BOTTOM QUARTILE OF PATIENTS BASED ON FEATURE SCORE CALCULATION
        # FROM INPUT GENE LIST
        ui.nav_panel("Gene-Signature Score Ranked DGE",
            ui.layout_sidebar(
                ui.sidebar(
                    ui.h4("Gene-Signature Score Ranked DGE"),
                    ui.accordion(

                        ui.accordion_panel("Genes for signature score",
                            ui.input_selectize("fs_dge_genes",
                                            "",
                                            [],
                                            multiple=True
                            ),
                        ),
                        ui.accordion_panel("Genes to highlight",
                            ui.input_selectize("fs_dge_highlight_genes",
                                            "",
                                            [],
                                            multiple=True
                            ),
                        ),
                        ui.accordion_panel("Compartment for gene-signature score calculation",
                            ui.input_selectize("fs_dge_compartment",
                                            "",
                                            qsh.compartments
                            ),
                        ),
                        ui.accordion_panel("Quartile",
                            ui.input_slider("fs_dge_slider",
                                            "Selects the top and bottom Quartile% of patients",
                                            value = .2,
                                            min = 0.05, max = .5,
                                            step = .05
                            ),
                        ),
                        ui.accordion_panel("DGE Compartment",
                            ui.input_selectize("fs_dge_compartment_for_dge",
                                            "DEGs calculated within this compartment",
                                            qsh.compartments
                            ),
                        ),
                        ui.accordion_panel("Thresholds",
                            ui.input_numeric("fs_dge_fc_thresh",
                                            "Fold Change Magnitude",
                                            value = 2
                            ),
                            ui.input_numeric("fs_dge_p_thresh",
                                            "-Log10(P-Value)",
                                            value=1),
                        )
                    ),
                    ui.input_action_button("fs_dge_run","RUN",icon=icon_svg("arrow-right")),
                    bg=panel_color
                ),
                ui.card(
                    ui.card(output_widget("gfs_ranked_dge")),
                    ui.layout_column_wrap(
                        ui.card(ui.card_header("Negative DEGs"),
                                ui.card_body(ui.output_data_frame("gfs_ranked_dge_bot")),
                                max_height=500),
                        ui.card(ui.card_header("Positive DEGs"),
                                ui.card_body(ui.output_data_frame("gfs_ranked_dge_top")),
                                max_height=500),
                    ),
                    ui.download_button("download_signaturecore_dges", "Download DEGs", icon=icon_svg("download")),
                ),
                ui.card(
                    ui.card_header("Expression levels for genes of interest"),
                    ui.card_body(output_widget('gfs_ranked_dge_highlighted_exprn_box'))
                ),
                bg=panel_color
            ),

            # METHODS
            ui.h4("Differential Gene Expression"),
            ui.p("""The Wilcoxon Rank-sum test is used to calculate differentially expressed genes.
                    P-values are corrected using Benjamini/Hochberg multiple testing correction.
                    """),
            ui.h5("Gene-Signature Score Calculation"),
            ui.p("""
                Data is subset to include archetyped patients only within a given compartment.
                Gene-signature scores are calculated by first calculating a per-gene Z-Score across patients
                and subsequently averaging the Z-Scores within each patient resulting in each patient being
                assigned a gene-signature score depending on their expression of the input gene set.

                """)
        ),
    ), # END DGE MENU

    ui.nav_control(ui.a("GitHub", href="https://github.com/hwismer/QuIPI", class_="nav-link",target="_blank")), # GITHUB LINK FOR ISSUES
    id = "quipi_top_nav",

    # CSS code controlling navbar display, background color etc.
    header=ui.tags.div(
        ui.tags.style(
             """
            .navbar-nav {
              display: flex;
              align-items: center; /* This is the key line for vertical centering */
              width: 100%;
              list-style: none;
              padding: 0;
            }
            .nav-item {
              margin: 0 15px;
            }
            body {
              background-color: #BDD0D5; /* Your desired navbar background color */
            }
            .nav-link { font-size: 16px;
                    color: black;
            }
            """
        ),
    ),
)

def server(input, output, session):

    quipi_tabs_visited = []

    ##### QUIPI

    @render_widget
    def cancer_glossary():
        return qsh.plot_cancer_glossary_table()
    
    @render_widget
    def plot_indication_breakdown():
        return qsh.plot_indication_breakdown()
    
    @render_widget
    def pancan_archetypes_home():
        fig = pp.plot_pancan_archetypes()
        fig.update_layout(template = "simple_white")
        return fig
    
    @render_widget
    def pancan_archetype_breakdown():
        return qsh.plot_archetype_beakdown()
    
    @render_widget
    @reactive.event(input.gene_factor_run)
    def pancan_archetypes_gfs():
        fig = pp.plot_pancan_archetypes()
        fig.update_layout(autosize=False, width=600, height=450,template = "simple_white")
        return fig
        
    @render_widget
    @reactive.event(input.pancan_umap_run)
    def pancan_subplots():
        
        transform = input.pancan_umap_transformation()
        genes = input.pancan_gene_input()
        compartment = input.pancan_compartment_input()

        fig =pp.plot_pancan_exprn_subplots(transform, genes, compartment)
        return fig


    ###### BOX / VIOL PLOTS

    # Reactive calculation where if multiple genes are selected, make the user choose a compartment for factor score calculation.
    @reactive.calc
    def box_viol_multi_options():
        if len(input.box_viol_gene_input()) > 1:
            choices = qsh.compartments
        else:
            choices = []
        return choices
        
    @render.ui
    def box_viol_multiple_genes():
        choices = box_viol_multi_options()
        if len(choices) > 1:
            return ui.input_selectize("box_viol_multiple_compartment", 
                                      "**Multiple genes selected. Choose compartment for factor score calculation.**", 
                                      choices,
                                      multiple=True)
        else:
            return None
        
    @reactive.effect
    @reactive.event(input.box_viol_x_category)  # Trigger when category changes
    def update_box_viol_selectize():
        x_cat = input.box_viol_x_category()
        new_options = qsh.categorials_opts_dict[x_cat]
        ui.update_selectize("box_viol_x_cat_filter", choices=new_options, selected=new_options)

    
    @render_widget
    @reactive.event(input.box_viol_run)
    def expression_box_viol():

        transform = input.box_viol_transformation()
        x_cat = qsh.categoricals_dict[input.box_viol_x_category()]
        x_cat_filts = input.box_viol_x_cat_filter()
        gene = input.box_viol_gene_input()
        group = input.box_viol_groupby()

        if len(gene) > 1:
            compartment_multiple = input.box_viol_multiple_compartment()
        else:
            compartment_multiple = None

        fig = bv.box_viol_exprn(transform, x_cat, x_cat_filts, gene, group, compartment_multiple)

        return fig
    
    ##### RIDGELINE
    #@render.plot() # height=100 * 10
    #@reactive.event(input.ridgeline_run)
    @reactive.Calc
    def ridgeline_helper():

        '''

        transform = input.box_viol_transformation()
        x_cat = qsh.categoricals_dict[input.box_viol_x_category()]
        x_cat_filts = input.box_viol_x_cat_filter()
        gene = input.box_viol_gene_input()
        group = input.box_viol_groupby()

        if len(gene) > 1:
            compartment_multiple = input.box_viol_multiple_compartment()
        else:
            compartment_multiple = None

        fig = bv.box_viol_exprn(transform, x_cat, x_cat_filts, gene, group, compartment_multiple)
        '''

        fig = bv.plot_ridgeline()
        print(fig)

        return fig, f"{1500}px"
    
        #height_param = max(750,15 * len(qsh.process_gene_text_input(input.gex_heatmap_text_genes())) )
        #height = f"{height_param}px"

        #return bv.plot_heatmap(genes, groupby, groups, transform), height

    @output(id="ridgeline_plot")
    @render.plot
    def ridgeline_plot():
        return ridgeline_helper()[0]

    @output(id="ridgeline")
    @render.ui
    @reactive.event(input.ridgeline_run)
    def run_ridgeline():
        height = ridgeline_helper()[1]
        return ui.output_plot("ridgeline_plot", height=height)
    




    ##### HEATMAP
    
    @reactive.effect
    @reactive.event(input.gex_heatmap_groupby)  # Trigger when category changes
    def update_heatmap_group_selectize():
        x_cat = input.gex_heatmap_groupby()
        cat_opts = qsh.categorical_choices[x_cat]
        ui.update_selectize("gex_heatmap_groupby_subset", choices=cat_opts, selected=cat_opts)


    @reactive.effect
    @reactive.event(input.gex_heatmap_splitby)  # Trigger when category changes
    def update_heatmap_group_selectize():
        split_cat = input.gex_heatmap_splitby()
        if split_cat != "---":
            cat_opts = qsh.categorical_choices[split_cat]
            ui.update_selectize("gex_heatmap_splitby_subset", choices=cat_opts, selected=cat_opts)
        else:
            ui.update_selectize("gex_heatmap_splitby_subset", choices=[])

    @reactive.Calc
    @reactive.event(input.gex_heatmap_run)
    def gex_heatmap_helper():
        genes = qsh.process_gene_text_input(input.gex_heatmap_text_genes())
        groupby = input.gex_heatmap_groupby()
        groups = input.gex_heatmap_groupby_subset()
        transform = input.gex_heatmap_transform()
        height_param = max(750,15 * len(qsh.process_gene_text_input(input.gex_heatmap_text_genes())) )
        height = f"{height_param}px"

        return bv.plot_heatmap(genes, groupby, groups, transform), height

    @output(id="gex_heatmap_plot")
    @render.plot
    def gex_heatmap_plot():
        return gex_heatmap_helper()[0]

    @output(id="gex_heatmap")
    @render.ui
    @reactive.event(input.gex_heatmap_run)
    def plot_container():
        height = gex_heatmap_helper()[1]
        return ui.output_plot("gex_heatmap_plot", height=height)

        
    
    ##### GENE EXPRESSION QUERY

    @reactive.calc
    @reactive.event(input.gene_expr_query_run)
    def gene_expr_query_backend():
        genes = qsh.process_gene_text_input(input.gene_expr_query_genes())
        indications = input.query_indication()
        compartments = input.query_compartment()
        archetypes = input.query_archetype()
        transform = input.query_transform()

        if transform == "TPM":
            df = pd.read_feather("./quipi_data/quipi_raw_tpm.feather", columns=qsh.non_genes + list(genes))
        elif transform == "Log2(TPM)":
            df = pd.read_feather("./quipi_data/quipi_log2_tpm.feather", columns=qsh.non_genes + list(genes))

        subset = df[(df["indication"].isin(indications)) & (df["compartment"].isin(compartments)) & (df["archetype"].isin(archetypes))]
        subset = subset[["patient","indication","sample_type","compartment","archetype"] + list(genes)]

        return subset
    
    @render.data_frame
    def gene_expr_query():
        return gene_expr_query_backend()
    
    @render.download(filename="quipi_query.csv")
    def download_query_table():
        df = gene_expr_query_backend()
        csv_buffer = StringIO()
        df.to_csv(csv_buffer, index=False)
        yield csv_buffer.getvalue()


    ##########  CORRELATION PLOTS

    ##### CORRELATION HEATMAP
    @render_widget
    @reactive.event(input.corr_run)
    def gene_correlation_heatmap():

        genes = input.corr_gene_input()
        indications = input.corr_indication()
        method = qsh.corr_methods[input.corr_method_input()]
        compartments = input.corr_compartment()
        archetypes = input.corr_archetype()
        tissues = [qsh.tissue_dict[tis] for tis in input.corr_tissue()]
        transform = input.corr_transform()

        return corr.gene_correlation_heatmap(genes, indications, method, compartments, archetypes, tissues, transform)
    

    ##### CROSS-COMPARTMENT CORRELATION HEATMAP
    @render_widget
    @reactive.event(input.comp_corr_mat_run)
    def comp_corr_mat():
        genes = list(input.comp_corr_mat_genes())
        comp1 = input.comp_corr_mat_compartment1()
        comp2 = input.comp_corr_mat_compartment2()
        transform = input.comp_corr_mat_transform()
        method = input.comp_corr_mat_method()

        indications = input.comp_corr_mat_indication()
        tissues = [qsh.tissue_dict[tis] for tis in input.comp_corr_mat_tissue()]
        archetypes = input.comp_corr_mat_archetype()

        fig = corr.compartment_correlation_heatmap(genes, comp1,comp2,transform,method,indications,tissues,archetypes)

        return fig


    ##### CATEGORICAL CORRELATION

    @render.data_frame
    @reactive.event(input.corr_cat_run)  
    async def corr_cat_gene_correlations(): 
        
        with ui.Progress(min=1, max = len(qsh.quipi_all_columns)) as p:
            p.set(message="Calculating", detail="Please Wait")
            genes = input.corr_cat_gene_input()
            category = input.corr_cat_category_input()
            categories = input.corr_cat_category_opts()
            range = input.corr_cat_slider()

            df = corr.categorical_correlation_table(genes, category, categories, range,p)
                
            return df
        

    ##### CROSS COMPARTMENT CATEGORICAL

    @render.data_frame
    @reactive.event(input.cross_comp_corr_cat_run)
    async def cross_comp_corr_cat_gene_correlations():
        genes = input.cross_comp_corr_cat_gene_input()
        comp1 = input.cross_comp_corr_comp1_input()
        comp2 = input.cross_comp_corr_comp2_input()
        range_slider = input.cross_comp_corr_cat_slider()
        method = input.cross_comp_corr_mat_method()
        transform = input.cross_comp_corr_mat_transform()

        with ui.Progress(min=0, max = len(qsh.quipi_all_columns) * len(genes)) as p:
            p.set(message="Calculating - I'm Accurate!", detail="This could take a while.")

            df = corr.cross_compartment_correlation_table(genes, comp1, comp2, range_slider, transform, method, p)

            return df

    @reactive.effect
    def corr_cat_update_choices():
        choices = input.corr_cat_category_input()
        ui.update_selectize("corr_cat_category_opts", choices=qsh.categorical_choices[choices])
    

    ##### PANCAN GENE-SIGNATURE OVERLAY

    @render_widget
    @reactive.event(input.gene_factor_run)
    def gene_factor_analysis():
        
        gene_set = list(input.gene_factor_genes())
        compartment = input.gene_factor_compartment()
        log2_subset_full = gf.calculate_gene_factor_score(gene_set, compartment)

        fig = px.scatter(log2_subset_full, x = "x_umap1", y = "x_umap2", 
                         color = "factor_score", color_continuous_scale="viridis",
                         labels = {"factor_score" : "Factor Score"})

        fig.update_layout(template="simple_white",autosize=False, width=600, height=450,legend_title_text = "Gene-Signature Score")
        
        fig.update_traces(marker=dict(size=12))
        fig.update_layout(legend_title_text = "Archetype")
        fig.update_yaxes(visible=False)
        fig.update_xaxes(visible=False)
        return fig
    


    ##### DGE BY COMPARTMENT SCORE

    @reactive.calc
    @reactive.event(input.dge_run)
    def feature_ranked_dge_reactive():

        fig, sig_pos, sig_neg, highlighted_boxplot = rpd.feature_ranked_dge(input.flow_score_to_rank(), 
                                      input.dge_compartment(), 
                                      input.dge_slider(),
                                      input.dge_fc_thresh(),
                                      input.dge_p_thresh(),
                                      list(input.dge_highlight_genes()))

        return fig, sig_pos, sig_neg, highlighted_boxplot

    
    @render_widget
    def compartment_featurescore_dge():
        return feature_ranked_dge_reactive()[0]
    
    @render.data_frame
    def compartment_featurescore_dge_top():
        return feature_ranked_dge_reactive()[1]
    
    @render.data_frame
    def compartment_featurescore_dge_bot():
        return feature_ranked_dge_reactive()[2]
    
    @render_widget
    def compartment_featurescore_dge_boxplot():
        boxplot = feature_ranked_dge_reactive()[3]
        if boxplot is not None:
            return boxplot
    
    @render.download(filename="FeatureScore_DEGs.csv")
    def download_featurescore_dges():
        df = pd.concat([feature_ranked_dge_reactive()[1],feature_ranked_dge_reactive()[2]], axis = 0, ignore_index=True)
        csv_buffer = StringIO()
        df.to_csv(csv_buffer, index=False)
        yield csv_buffer.getvalue()



    ##### DGE BY GENE-SIGNATURE
    @reactive.calc
    @reactive.event(input.fs_dge_run)
    def factor_ranked_dge_reactive():
        
        if len(input.fs_dge_genes()) != 0:
            fig, sig_pos, sig_neg, highlight_boxplot = rpd.factor_ranked_dge(list(input.fs_dge_genes()),
                                                        input.fs_dge_compartment(),
                                                        input.fs_dge_slider(),
                                                        input.fs_dge_compartment_for_dge(),
                                                        input.fs_dge_fc_thresh(),
                                                        input.fs_dge_p_thresh(),
                                                        list(input.fs_dge_highlight_genes()))

            return fig, sig_pos, sig_neg, highlight_boxplot
    
    @render_widget
    def gfs_ranked_dge():
        return factor_ranked_dge_reactive()[0]

    
    @render.data_frame
    def gfs_ranked_dge_top():
        return factor_ranked_dge_reactive()[1]
    
    @render.data_frame
    def gfs_ranked_dge_bot():
        return factor_ranked_dge_reactive()[2]
    
    @render_widget
    def gfs_ranked_dge_highlighted_exprn_box():
        boxplot = factor_ranked_dge_reactive()[3]
        
        if boxplot != None :
            return boxplot
    
    @render.download(filename="SignatureScore_DEGS.csv")
    def download_signaturecore_dges():
        df = pd.concat([factor_ranked_dge_reactive()[1],factor_ranked_dge_reactive()[2]], axis = 0, ignore_index=True)
        csv_buffer = StringIO()
        df.to_csv(csv_buffer, index=False)
        yield csv_buffer.getvalue()

    ## HELPER: Populate individual gene selection boxes to avoid long startup.
    @reactive.effect
    def populate_gene_selections():
        current_tab = input.quipi_top_nav()
        if current_tab not in quipi_tabs_visited:
            quipi_tabs_visited.append(current_tab)
            if current_tab in quipi_tabs_mapped_to_gene_inputs:
                for id in quipi_tabs_mapped_to_gene_inputs[current_tab]:
                    ui.update_selectize(
                        id,
                        choices=qsh.genes,
                        selected=[],
                        server=True,
                    )


# Create the Shiny app
app_dir = Path(__file__).parent
app = App(app_ui, server,static_assets= app_dir / "www")



