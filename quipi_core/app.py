from shiny import App, render, ui, reactive
from shiny.types import ImgData
from shinyswatch import theme
import plotly.express as px
from shinywidgets import output_widget, render_widget 

import quipi_shared as qsh
import humu_shared as hsh

import numpy as np
import pandas as pd
from io import StringIO

import ranked_patient_dge as rpd
import gene_factor as gf
import box_viol_expression_plot as bv
import pancan_plots as pp
import correlation_analysis as corr
from pathlib import Path

import humu_flow_boxplot as hfb
import humu_gex_plots as hxp


quipi_tabs_mapped_to_gene_inputs = {"Box/Violin Plots" : ["box_viol_gene_input"],
                              "Gene Expresssion Bar Plots" : ["gene_expr_bar_genes"],
                              "Dotplots" : ["gex_dot_genes"],
                              "Correlation Matrix" : ["corr_gene_input"],
                              "Compartment Correlation Matrix" : ["comp_corr_mat_genes"],
                              "One-Vs-All Correlation Table" : ["corr_cat_gene_input"],
                              "PanCan UMAP Gene Expression" : ["pancan_gene_input"],
                              "PanCan Gene-Signature Overlay" : ["gene_factor_genes"],
                              "Feature Score Ranked DGE" : ["dge_highlight_genes"],
                              "Gene-Signature Score Ranked DGE" : ["fs_dge_genes","fs_dge_highlight_genes"],
                              "Query Gene Expression" : ["gene_expr_query_genes"],
                              "Cross-Compartment Correlation Table": ["cross_comp_corr_cat_gene_input"],

}

humu_tabs_mapped_to_gene_inputs = {"Mouse Gene Expression Box Plots" : [["humu_gex_box_gene", "Mouse"]],
                                   "Mouse Gene Expression Dotplots" : [["humu_gex_dot_gene", "Mouse"]],
                                   "HuMu Gene Expression Comparison" : [["humu_box_comp_human_genes", "Human"],["humu_box_comp_mu_genes", "Mouse"]],
}

panel_color = "#f0f0f0"

# Define the UI
app_ui = ui.page_fluid(
    ui.tags.style("""
            body { background-color: #b8d1d6; }  /* Background color */
            .nav-link { font-size: 20px;
                        color: white;
            }
            .navbar-nav .nav-link.active {
                color: white !important;  /* Keep text white when selected */
            }
            .my-image {
                border: 3px solid #333; /* Adjust border width, style, and color as needed */
            }
    """),

    ui.navset_bar(

        ui.nav_panel("QuIPI Human",
            ui.tags.div(
                ui.tags.div(
                    ui.tags.img(src="quipi.png", style="height: 100px; margin-left: 10px; padding-top: 10px"),  # Left-aligned image
                    ui.tags.span("", style="font-size: 50px; font-weight: bold;"),
                    style="display: flex; align-items: flex-end; gap: 10px;" # Flexbox for horizontal alignment
                ),
                style="padding-bottom: 5px;"
            ),

            ui.page_navbar(
            
                ## HOME PAGE
                ui.nav_panel("Home",
                    ui.panel_title("QuIPI - Querying IPI"),
                    ui.p("Here are some useful references."),
                    ui.layout_column_wrap(
                        ui.card(ui.card_header("Cancer Indication Breakdown"),
                                output_widget("plot_indication_breakdown"), 
                                full_screen=False,
                                ),
                        ui.card(ui.card_header("Cancer Indication Nomenclature"),
                            output_widget("cancer_glossary"),
                            full_screen=False),
                        ui.card(ui.card_header("PanCan Archetype UMAP"),
                                ui.card_body(output_widget("pancan_archetypes_home")),
                                ui.card_footer(
                                    ui.HTML(f'<a href="{"https://pubmed.ncbi.nlm.nih.gov/34963056/"}" target="_blank">Discovering dominant tumor immune archetypes in a pan-cancer census. Combes AJ, Samad B, Tsui J, et al. Cell. 2022</a>')
                                ),
                                full_screen=False,     
                        ),
                        ui.card(ui.card_header("PanCan Archetype Breakdown"),
                                ui.card_body(output_widget("pancan_archetype_breakdown")),
                                ui.card_footer("""
                                            Nota Bene: Unclassified patients do not have Feature Scores assigned to them.
                                            Any calculations involving a feature score are subset to include archetyped patients only.
                                            """)),
                        width=.5
                    ),
                ),

                ##### GENE EXPRESSION
                ui.nav_menu("Gene Expression",
                    
                    ui.nav_panel("Box/Violin Plots",
                                 
                        ui.layout_sidebar(

                            ui.sidebar(

                                ui.h4("Explore gene expression by category"),
                                ui.input_action_button("box_viol_run", "Run"),
                                ui.input_selectize("box_viol_plot",
                                                    "Select Plot Type:",
                                                    ["Boxplot","Violin Plot"]),
                                ui.input_selectize("box_viol_gene_input",
                                                    "Select Genes:",
                                                    [],
                                                    multiple=True),
                                ui.output_ui("box_viol_multiple_genes"),
                                ui.input_selectize("box_viol_x_category",
                                                    "Select X-Axis Category:",
                                                    list(qsh.categoricals_dict.keys()),
                                                    selected = "Compartment"),
                                ui.input_selectize("box_viol_x_cat_filter",
                                                "**Subset X-Axis Categories.**",
                                                [],
                                                multiple=True),
                                ui.input_selectize("box_viol_groupby",
                                                    "Group by:",
                                                    ["---"] + list(qsh.categoricals_dict.keys()),
                                                    selected="---"),
                                ui.input_selectize("box_viol_transformation",
                                                    "Select TPM Transformation: ",
                                                    ["TPM", "Log2(TPM)"],
                                                    multiple= False,
                                                    selected= "Log2(TPM)"),
                                bg=panel_color
                            ),
                        
                        ui.card(ui.card_body(output_widget("expression_box_viol")),
                                ui.card_footer("Click button in the bottom right for fullscreen view."),
                                full_screen=True),
                        bg=panel_color
                        ),
                    ),

                    ui.nav_panel("Dotplots",
                        ui.layout_sidebar(
                            ui.sidebar(
                                ui.input_selectize("gex_dot_genes", "Choose Genes to plot:", [], multiple=True),
                                ui.card(
                                    ui.input_selectize("gex_dot_groupby", "Group by:", list(qsh.categorical_choices.keys())),
                                    ui.input_selectize("gex_dot_groupby_subset", "Subset Groupby Categories:", [], multiple=True),
                                ),
                                ui.card(
                                    ui.input_selectize("gex_dot_splitby", "Split by:", ["---"] + list(qsh.categorical_choices.keys()), selected="---"),
                                    ui.input_selectize("gex_dot_splitby_subset", "Subset Splitby Categories:", [], multiple=True), 
                                ),
                                ui.input_selectize("gex_dot_transform", "Select TPM Transformation:", ["Log2(TPM)", "TPM"], selected = "TPM"),
                                ui.input_switch("gex_dot_swap", "Swap Axes"),
                                ui.input_action_button("gex_dot_run", "RUN"),
                                bg=panel_color
                            ),
                            ui.card(ui.output_plot("gex_dotplot", fill=True), full_screen=True)   
                        )
                    ),

                    ui.nav_panel("Query Gene Expression",
                        ui.layout_sidebar(
                            ui.sidebar(
                                ui.h4("Subset and download gene expression data"),
                                ui.input_action_button("gene_expr_query_run", "Run"),
                                ui.download_button("download_query_table", "Download CSV"),
                                ui.input_selectize("gene_expr_query_genes",
                                                "Select Genes:",
                                                [],
                                                multiple=True,
                                                options = {"server":True}),
                                ui.input_selectize("query_indication",
                                                    "Select Indications:",
                                                    choices=qsh.indications,
                                                    selected=qsh.indications,
                                                    multiple = True),
                                ui.input_selectize("query_compartment",
                                                    "Select Compartments:",
                                                    choices=qsh.compartments,
                                                    selected=qsh.compartments,
                                                    multiple=True),
                                ui.input_selectize("query_archetype",
                                                    "Select Archetypes:",
                                                    choices=qsh.archetypes,
                                                    multiple=True,
                                                    selected=qsh.archetypes),
                                ui.input_selectize("query_transform",
                                                    "Select TPM Transformation:",
                                                    choices=["TPM", "Log2(TPM)"],
                                                    selected="Log2(TPM)"),
                                bg=panel_color
                            ),
                            ui.card(ui.output_data_frame("gene_expr_query")),
                            #ui.download_button("download_query_table", "Download CSV"),
                            bg=panel_color
                        ),
                    )
                ), # END GENE EXPRESSION MENU

                # Pairwise correlation plot between user-defined genes.
                ui.nav_menu("Correlation Plots",
                            
                    ui.nav_panel("Correlation Matrix",
                        ui.h4("Explore the correlation of genes across the entire IPI dataset."),
                        ui.layout_sidebar(
                            ui.sidebar(
                                ui.input_action_button("corr_run", "Run"),
                                ui.input_selectize("corr_gene_input",
                                                    "Select Genes:",
                                                    [],
                                                    multiple=True,
                                                    options = {'server':True}),
                                ui.input_selectize("corr_indication",
                                                    "Select Indications:",
                                                    choices=qsh.indications,
                                                    selected=qsh.indications,
                                                    multiple = True),
                                ui.input_selectize("corr_tissue",
                                                    "Select Tissue:",
                                                    choices=["Tumor", "Normal"],
                                                    selected = ["Tumor", "Normal"],
                                                    multiple = True),
                                ui.input_selectize("corr_compartment",
                                                    "Select Compartments:",
                                                    choices=qsh.compartments,
                                                    selected=qsh.compartments,
                                                    multiple=True),
                                ui.input_selectize("corr_archetype",
                                                    "Select Archetypes:",
                                                    choices=qsh.archetypes,
                                                    multiple=True,
                                                    selected=qsh.archetypes),
                                ui.input_selectize("corr_transform",
                                                    "Select TPM Transformation:",
                                                    choices=["TPM", "Log2(TPM)"],
                                                    selected="Log2(TPM)"),
                                ui.input_selectize("corr_method_input",
                                                    "Select Method:",
                                                    choices=["Pearson", "Spearman"],
                                                    selected="Spearman"),
                                bg=panel_color
                            ),
                            ui.card(ui.card_body(output_widget("gene_correlation_heatmap")),
                                    ui.card_footer("Click button in the bottom right for fullscreen view."),
                                    full_screen=True),
                            bg=panel_color
                        )
                    ),

                    ui.nav_panel("Compartment Correlation Matrix",
                        ui.h4("Explore the correlation of chosen genes across two compartments."),
                        ui.layout_sidebar(
                            ui.sidebar(
                                ui.input_action_button("comp_corr_mat_run", "Run"),
                                ui.input_selectize("comp_corr_mat_genes",
                                                "Select Genes:",
                                                [],
                                                multiple=True,
                                                options = {"server":True}),
                                ui.input_selectize("comp_corr_mat_compartment1",
                                                    "Select Compartment A:",
                                                    choices=qsh.compartments),
                                ui.input_selectize("comp_corr_mat_compartment2",
                                                "Select Compartment B:",
                                                choices=qsh.compartments),
                                ui.input_selectize("comp_corr_mat_method",
                                                "Select Correlation Method:",
                                                choices = ["Pearson", "Spearman"],
                                                selected="Spearman"),
                                ui.input_selectize("comp_corr_mat_transform",
                                                    "Select TPM Transformation:",
                                                    choices=["TPM", "Log2(TPM)"],
                                                    selected="Log2(TPM)"),
                                ui.input_selectize("comp_corr_mat_indication",
                                                "Select Indications:",
                                                choices=qsh.indications,
                                                selected=qsh.indications,
                                                multiple = True),
                                ui.input_selectize("comp_corr_mat_tissue",
                                                    "Select Tissue:",
                                                    choices=["Tumor", "Normal"],
                                                    selected = ["Tumor", "Normal"],
                                                    multiple = True),
                                ui.input_selectize("comp_corr_mat_archetype",
                                                    "Select Archetypes:",
                                                    choices=qsh.archetypes,
                                                    multiple=True,
                                                    selected=qsh.archetypes),
                                bg=panel_color
                            ),
                            ui.card(ui.card_body(output_widget("comp_corr_mat")),
                                    ui.card_footer("Click button in the bottom right for fullscreen view."),
                                    full_screen=True),
                            bg=panel_color
                        )
                    ),

                    ui.nav_panel("One-Vs-All Correlation Table",
                        ui.h3("One-Vs-All correlation analysis within a chosen category."),
                        ui.layout_sidebar(
                                ui.sidebar(
                                    ui.input_selectize("corr_cat_gene_input",
                                                        "Select Gene:",
                                                        [],
                                                        multiple=False,
                                                        options = {'server':True}),
                                    ui.input_selectize("corr_cat_category_input",
                                                    "Select Category:",
                                                    choices=list(qsh.categorical_choices.keys())),
                                    ui.input_selectize("corr_cat_category_opts",
                                                    "Select Subcategories:",
                                                    choices = [],
                                                    multiple=True),
                                    ui.input_slider("corr_cat_slider", "Select Correlation Coefficient Range", min=-1,max=1, value = (0,1), step=.1),
                                    ui.input_action_button("corr_cat_run", "Run"),
                                    bg=panel_color           
                                ),
                            ui.output_data_frame("corr_cat_gene_correlations"),
                            bg=panel_color
                        )
                    ),
                    
                    ui.nav_panel("Cross-Compartment Correlation Table",
                        ui.h3("Cross-Compartment correlation for a chosen gene."),
                        ui.layout_sidebar(
                                ui.sidebar(
                                    ui.input_selectize("cross_comp_corr_cat_gene_input",
                                                        "Select Gene:",
                                                        [],
                                                        multiple=True,
                                                        options = {'server':True}),
                                    ui.input_selectize("cross_comp_corr_comp1_input",
                                                    "Select Compartment 1:",
                                                    choices=list(qsh.compartments)),
                                    ui.input_selectize("cross_comp_corr_comp2_input",
                                                    "Select Compartment 2",
                                                    choices = list(qsh.compartments),
                                                    multiple=True),
                                    ui.input_slider("cross_comp_corr_cat_slider", "Select Correlation Coefficient Range", min=-1,max=1, value = (0,1), step=.1),
                                    ui.input_selectize("cross_comp_corr_mat_method",
                                                "Select Correlation Method:",
                                                choices = ["Pearson", "Spearman"],
                                                selected="Spearman"),
                                    ui.input_selectize("cross_comp_corr_mat_transform",
                                                    "Select TPM Transformation:",
                                                    choices=["TPM", "Log2(TPM)"],
                                                    selected="Log2(TPM)"),
                                    ui.input_action_button("cross_comp_corr_cat_run", "Run"),
                                    bg=panel_color
                                ),
                            ui.output_data_frame("cross_comp_corr_cat_gene_correlations"),
                            bg=panel_color
                        )
                    )
                ), # END CORRELATION MENU

                ui.nav_menu("PanCan Archetype UMAP",
                    # Tab where user can select multiple eachs and view their expression overlayed on the individual
                    # PanCan UMAPs
                    ui.nav_panel("PanCan UMAP Gene Expression",
                        ui.layout_sidebar(
                            ui.sidebar(
                                ui.h4("PanCan UMAP Gene Expression"),
                                ui.input_action_button("pancan_umap_run", "Run"),
                                ui.input_selectize("pancan_gene_input",
                                                    "Select Genes:",
                                                    [],
                                                    multiple = True),
                                ui.input_selectize("pancan_compartment_input",
                                                "Select Compartment:",
                                                qsh.compartments
                                                ),
                                ui.input_selectize("pancan_umap_transformation",
                                                "Select TPM Transformation:",
                                                ["TPM", "Log2(TPM)"],
                                                multiple= False,
                                                selected= "Log2(TPM)"),
                                bg=panel_color
                            ),
                            output_widget("pancan_subplots"),
                            bg=panel_color
                        ),
                    ),

                    # Given user-defined gene list, plots the gene factor score on top of the PanCan UMAP
                    ui.nav_panel("PanCan Gene-Signature Overlay",
                        ui.layout_sidebar(
                            ui.sidebar(
                                ui.h4("PanCan Gene-Signature Overlay"),
                                ui.input_action_button("gene_factor_run", "Run"),
                                ui.input_selectize("gene_factor_genes",
                                                    "Select Genes:",
                                                    [],
                                                    multiple=True),
                                ui.input_selectize("gene_factor_compartment",
                                                "Select Compartment:",
                                                qsh.compartments),
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
                ui.nav_menu("Ranked DGE",
                            
                    ui.nav_panel("Feature Score Ranked DGE",
                        ui.layout_sidebar(
                            ui.sidebar(
                                ui.h4("Feature Score Ranked DGE"),
                                ui.input_action_button("dge_run", "Run"),
                                ui.input_selectize("flow_score_to_rank",
                                                "Feature Score For Ranking:",
                                                list(qsh.feature_scores.keys())),
                                ui.input_selectize("dge_highlight_genes",
                                                "Genes of interest to be highlighted.",
                                                [],
                                                    multiple=True),
                                ui.input_slider("dge_slider",
                                                "Quartile:",
                                                value = .2,
                                                min = 0.05,
                                                max = .5,
                                                step = .05),
                                ui.input_selectize("dge_compartment",
                                                    "DGE Compartment:",
                                                    qsh.compartments),
                                ui.input_numeric("dge_fc_thresh",
                                                "Log2FC Threshold:",
                                                value=2),
                                ui.input_numeric("dge_p_thresh",
                                                "-Log10(P-Value) Threshold:",
                                                value = 1),
                                bg=panel_color
                            ),
                            ui.card(output_widget("compartment_featurescore_dge")),
                            ui.layout_column_wrap(
                                ui.card(
                                    ui.card_header("Negative DGEs"),
                                    ui.card_body(ui.output_data_frame("compartment_featurescore_dge_bot")),
                                ),
                                ui.card(
                                    ui.card_header("Positive DGEs"),
                                    ui.card_body(ui.output_data_frame("compartment_featurescore_dge_top"))
                                ),
                            ),
                            ui.download_button("download_featurescore_dges", "Download DGEs",),
                            ui.card(
                                ui.card_header("Expression Levels for genes of interest"),
                                ui.card_body(output_widget("compartment_featurescore_dge_boxplot"))
                            ),
                            bg=panel_color
                        ),
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

                    # DGE between the top and bottom quantiles of patients based on
                    # gene factor scores calculated from a collection of genes.
                    ui.nav_panel("Gene-Signature Score Ranked DGE",
                        ui.layout_sidebar(
                            ui.sidebar(
                                ui.h4("Gene-Signature Score Ranked DGE"),
                                ui.input_action_button("fs_dge_run","Run"),
                                ui.input_selectize("fs_dge_genes",
                                                "Genes for gene-signature score:",
                                                [],
                                                multiple=True),
                                ui.input_selectize("fs_dge_highlight_genes",
                                                "Genes of interest to be highlighted.",
                                                [],
                                                    multiple=True),
                                ui.input_selectize("fs_dge_compartment",
                                                "Compartment to calculate gene-signature score:",
                                                qsh.compartments),
                                ui.input_slider("fs_dge_slider",
                                                "Quartile:",
                                                value = .2,
                                                min = 0.05, max = .5,
                                                step = .05),
                                ui.input_selectize("fs_dge_compartment_for_dge",
                                                "Compartment for DGE:",
                                                qsh.compartments),
                                ui.input_numeric("fs_dge_fc_thresh",
                                                "Fold Change Magnitude Threshold",
                                                value = 2),
                                ui.input_numeric("fs_dge_p_thresh",
                                                "-Log10(P-Value) Threshold",
                                                value=1),
                                bg=panel_color
                            ),
                            ui.card(output_widget("gfs_ranked_dge")),
                            ui.layout_column_wrap(
                                ui.card(ui.card_header("Negative DGEs"),
                                        ui.card_body(ui.output_data_frame("gfs_ranked_dge_bot")),
                                        max_height=500),
                                ui.card(ui.card_header("Positive DGEs"),
                                        ui.card_body(ui.output_data_frame("gfs_ranked_dge_top")),
                                        max_height=500),
                            ),
                            ui.download_button("download_signaturecore_dges", "Download DGEs"),
                            ui.card(
                                ui.card_header("Expression levels for genes of interest"),
                                ui.card_body(output_widget('gfs_ranked_dge_highlighted_exprn_box'))
                            ),
                            bg=panel_color
                        ),
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
                ), # END DGE MENU
            ),
                id = "quipi_top_nav",
                theme=theme.cosmo,
                bg = "#1a1807"
            ),
        ),



        ##### HUMU 


        ui.nav_panel("QuIPI HuMu",
            ui.tags.div(
                ui.tags.img(src="humu.png", style="width: 90px; margin-left: 10px; margin-top: 21px; margin-bottom: 10px",class_="my-image"),  # Left-aligned image
                ui.tags.span("QuIPI HuMu", style="font-size: 40px; font-weight: bold;"),
                style="display: flex; align-items: flex-end; gap: 10px;" # Flexbox for horizontal alignment
            ),
            ui.page_navbar(
                ui.nav_panel("Home",
                    ui.layout_column_wrap(
                        ui.card(
                            ui.h4("The Human-to-Mouse Cancer Translator Project (HuMu) aims to immuno-profile a series of common and exceptional models of cancer in mice to benchmark them against the diversity of TMEs in Human cancer (i.e immune archetypes described in Combes, Samad, et al. Cell 2022). The HuMu dataset contains high-dimensional cytometry data (CyTOF) of 15 murine models that we use to study high level tumor-immune composition, as well as single-cell sequencing data from 9 of these models that we use to dissect more granular gene expression profiles across populations and tumor models"),
                        ),
                        ui.card(
                            ui.tags.img(src="quipi_humu_reference.png",style="width: 750px; margin-left: 10px; margin-top: 21px; margin-bottom: 10px"),  # Left-aligned image
                        )
                    )
                ),

                ui.nav_panel("Flow Boxplots",
                    #ui.h4("Explore gene expression by category"),
                    ui.layout_sidebar(
                        ui.sidebar(
                            ui.input_selectize("humu_box_score_1", "Choose First Score", hsh.flow_scores),
                            ui.input_selectize("humu_box_score_2", "Choose Second Score if ratio desired", ["---"] + hsh.flow_scores, selected="---"),
                            ui.input_selectize("humu_box_x_cat", "Select X-Axis Category", ["Species", "Group"]),
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

                ui.nav_panel("HuMu Gene Expression Comparison",
                    ui.layout_sidebar(
                        ui.sidebar(
                            ui.input_selectize("humu_box_comp_human_genes", "Choose Human Gene", []),
                            ui.input_selectize("humu_box_comp_mu_genes", "Choose Murine Gene", []),
                            ui.input_action_button("humu_box_comp_run", "RUN"),
                        bg=panel_color
                        ),
                        ui.card(output_widget("humu_gene_comparision")),
                        bg=panel_color
                    )
                ),

                ui.nav_panel("Mouse Gene Expression Box Plots",
                    ui.layout_sidebar(
                        ui.sidebar(
                            ui.input_selectize("humu_gex_box_gene", "Choose Gene to plot:", []),
                            ui.input_selectize("humu_gex_box_x_cat", "Choose X-Axis Category:", hsh.categoricals_opts),
                            ui.input_selectize("humu_gex_box_cat_subset", "Subset Categories:", [], multiple=True),
                            ui.input_selectize("humu_gex_box_groupby", "Group by:", ["---"] + hsh.categoricals_opts, selected="---"),
                            ui.input_selectize("humu_gex_box_splitby", "Split by:", ["---"] + hsh.categoricals_opts, selected="---"),
                            ui.input_action_button("humu_gene_box_run", "RUN"),
                            bg=panel_color
                        ),
                        ui.card(output_widget("gex_box"), full_screen=True),
                        bg=panel_color
                    )       
                ),

                ui.nav_panel("Mouse Gene Expression Dotplots",
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
                    ui.card(ui.output_plot("humu_plot_gex_dotplot"), min_height="750px", full_screen=True)
                    )       
                ),
            ui.nav_spacer(),
            id = "humu_top_nav",
            theme=theme.cosmo,
            bg = "#1a1807"
            )
        ),
    title = "",
    ),
    theme=theme.cosmo,
    bg = "#1a1807"
)

def server(input, output, session):

    quipi_tabs_visited = []
    humu_tabs_visited = []

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
        fig.update_layout(autosize=False, width=650, height=500,template = "simple_white")
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
        plot_type = input.box_viol_plot()

        if len(gene) > 1:
            compartment_multiple = input.box_viol_multiple_compartment()
        else:
            compartment_multiple = None

        fig = bv.box_viol_exprn(transform, x_cat, x_cat_filts, gene, group, plot_type, compartment_multiple)

        return fig
    
    ##### DOTPLOTS
    
    @reactive.effect
    @reactive.event(input.gex_dot_groupby)  # Trigger when category changes
    def update_dotplot_group_selectize():
        x_cat = input.gex_dot_groupby()
        cat_opts = qsh.categorical_choices[x_cat]
        ui.update_selectize("gex_dot_groupby_subset", choices=cat_opts, selected=cat_opts)


    @reactive.effect
    @reactive.event(input.gex_dot_splitby)  # Trigger when category changes
    def update_dotplot_group_selectize():
        split_cat = input.gex_dot_splitby()
        if split_cat != "---":
            cat_opts = qsh.categorical_choices[split_cat]
            ui.update_selectize("gex_dot_splitby_subset", choices=cat_opts, selected=cat_opts)
        else:
            ui.update_selectize("gex_dot_splitby_subset", choices=[])

    @render.plot
    @reactive.event(input.gex_dot_run)
    def gex_dotplot():
        genes = list(input.gex_dot_genes())
        groupby = input.gex_dot_groupby()
        groups = input.gex_dot_groupby_subset()
        splitby = input.gex_dot_splitby()
        splits = input.gex_dot_splitby_subset()
        transform = input.gex_dot_transform()
        swap = input.gex_dot_swap()

        return bv.dotplot(genes, groupby, groups, splitby, splits, transform, swap)

        
    
    ##### GENE EXPRESSION QUERY

    @reactive.calc
    @reactive.event(input.gene_expr_query_run)
    def gene_expr_query_backend():
        genes = input.gene_expr_query_genes()
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
    
    @render.download(filename="FeatureScore_DGES.csv")
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
    
    @render.download(filename="SignatureScore_DGES.csv")
    def download_signaturecore_dges():
        df = pd.concat([factor_ranked_dge_reactive()[1],factor_ranked_dge_reactive()[2]], axis = 0, ignore_index=True)
        csv_buffer = StringIO()
        df.to_csv(csv_buffer, index=False)
        yield csv_buffer.getvalue()





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
        fig = hfb.box_humu_flow(score1, score2, x_cat, x_cat_filter, group)
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
    def humu_gene_comparision():
        human_gene = input.humu_box_comp_human_genes()
        mouse_gene = input.humu_box_comp_mu_genes()

        fig = hxp.humu_box_comparison(human_gene, mouse_gene)

        return fig


    

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
                            choices=qsh.genes,
                            selected=[],
                            server=True,
                        )


# Create the Shiny app
app_dir = Path(__file__).parent
app = App(app_ui, server,static_assets= app_dir / "www")



