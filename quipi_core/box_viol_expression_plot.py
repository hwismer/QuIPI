import plotly.express as px
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

import shared as sh
import gene_factor as gf


def box_viol_exprn(transform, x_cat, x_cat_filts, genes, groupby, plot_type, compartment_multiple):


    if groupby in sh.categoricals_dict:
        group = sh.categoricals_dict[groupby]
    else:
        group = x_cat

    new_cats = dict(sh.categoricals_dict_reversed)
    new_cats.update(factor_score="Gene-Signature Score")

    color = sh.color_dict[group]

    if len(genes) > 1 and len(compartment_multiple) != 0:
        input_arr = gf.calculate_gene_factor_score_all_patients(list(genes), compartment_multiple)
        input_arr = input_arr[input_arr[x_cat].isin(x_cat_filts)]
        
        if plot_type == "Boxplot":
            fig = px.box(input_arr, x = x_cat, y = "factor_score", color = group,
                        color_discrete_map=color,
                        labels=new_cats)
        else:
            fig = px.violin(input_arr, x = x_cat, y ="factor_score", color = group,
                            color_discrete_map=color,
                            labels=new_cats)
        return fig

    elif len(genes) == 1:
        if transform == "TPM":
            input_arr = pd.read_feather("./quipi_data/quipi_raw_tpm.feather", columns=sh.non_genes + list(genes))
        elif transform == "Log2(TPM)":
            input_arr = pd.read_feather("./quipi_data/quipi_log2_tpm.feather", columns=sh.non_genes + list(genes))
    
        if len(genes) != 0:

            if plot_type == "Boxplot":
                fig = px.box(input_arr, x = x_cat,y = genes[0],color = group,
                            color_discrete_map=color,
                            labels=sh.categoricals_dict_reversed)
            else:
                fig = px.violin(input_arr, x = x_cat, y = genes[0], color = group,
                                color_discrete_map=color,
                                labels=sh.categoricals_dict_reversed)
                
            fig.update_layout(title_text= transform + genes[0], title_x=0.5)
            return fig
        
def dotplot(genes, groupby, groups, splitby, splits, transform, swap):
    adata = sc.read_h5ad("./quipi_data/quipi_adata.h5ad", backed="r")

    fig, ax = plt.subplots()
    
    groupby = sh.categoricals_dict[groupby]

    if transform == "TPM":
        transform = None
    else:
        transform = "log2"

    if splitby != "---":
        splitby = sh.categoricals_dict[splitby]
        adata = adata[adata.obs[splitby].isin(splits) & (adata.obs[groupby].isin(groups))]
    else:
        adata = adata[adata.obs[groupby].isin(groups)]

    if groupby == splitby:
        vars = groupby
    else:
        vars = [groupby, splitby] if splitby != "---" else groupby
    
    if swap:
        sc.pl.DotPlot(adata, var_names=genes, groupby=vars, ax=ax, layer=transform).swap_axes().make_figure()
    else:
        sc.pl.DotPlot(adata, var_names=genes, groupby=vars, ax=ax, layer=transform).make_figure()

    fig.tight_layout()
    plt.subplots_adjust(top=0.95, bottom=0.15, left=0.1, right=0.95, hspace=0.0, wspace=0.0)
    ax = plt.gca()
    return fig
