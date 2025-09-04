import plotly.express as px
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np

import quipi_shared as sh
import gene_factor as gf

import seaborn as sns


# BOXPLOT PLOTTING FUNCTION
def box_viol_exprn(transform, x_cat, x_cat_filts, genes, groupby, compartment_multiple):


    if groupby in sh.categoricals_dict:
        group = sh.categoricals_dict[groupby]
    else:
        group = x_cat

    orders = sh.quipi_orders[x_cat]

    new_cats = dict(sh.categoricals_dict_reversed)
    new_cats.update(factor_score="Gene-Signature Score")

    color = sh.color_dict[group]

    if len(genes) > 1 and len(compartment_multiple) != 0:

        input_arr = gf.calculate_gene_factor_score_all_patients(list(genes), compartment_multiple)
        input_arr = input_arr[input_arr[x_cat].isin(x_cat_filts)]
    
        fig = px.box(input_arr, x = x_cat, y = "factor_score", color = group,
                    color_discrete_map=color, category_orders=sh.quipi_orders,
                    labels=new_cats)
        return fig

    elif len(genes) == 1:

        if transform == "TPM":
            input_arr = pd.read_feather("./quipi_data/quipi_raw_tpm.feather", columns=sh.non_genes + list(genes))
        elif transform == "Log2(TPM)":
            input_arr = pd.read_feather("./quipi_data/quipi_log2_tpm.feather", columns=sh.non_genes + list(genes))

        input_arr = input_arr[input_arr[x_cat].isin(x_cat_filts)]

        if len(genes) != 0:

            fig = px.box(input_arr, x = x_cat,y = genes[0],color = group,
                        color_discrete_map=color, category_orders={x_cat: orders},
                        labels=sh.categoricals_dict_reversed)
                
            fig.update_layout(title_text= transform + genes[0], title_x=0.5)
            return fig
        

# DOTPLOT - CURRENTLY DEPRECATED REPLACED BY HEATMAP
def plot_dotplot(genes, groupby, groups, splitby, splits, transform, swap):

    CUTOFF=0

    if transform == "TPM":
        path = "./quipi_data/quipi_raw_tpm.feather"
    elif transform == "Log2(TPM)":
        path = "./quipi_data/quipi_log2_tpm.feather"

    if groupby == splitby:
        splitby = "---"

    groupby = sh.categoricals_dict[groupby]

    if splitby == "---":
        df = pd.read_feather(path, 
                         columns = genes + [groupby])
        df = df[df[groupby].isin(groups)]
        df["group_splits"] = df[groupby].astype(str)
        
    else:
        splitby = sh.categoricals_dict[splitby]
        df = pd.read_feather(path, 
                         columns = genes + [groupby] + [splitby])
        df = df[(df[groupby].isin(groups)) & (df[splitby].isin(splits))]
        df["group_splits"] = df[groupby].astype(str) + '_' + df[splitby].astype(str)

    
    
    expressed = df[genes] > CUTOFF
    expressed["group_splits"] = df["group_splits"]

    dot_size_sum = expressed.groupby("group_splits")[genes].sum()
    dot_size_total = expressed.groupby("group_splits")[genes].count()
    
    dot_size_df = dot_size_sum / dot_size_total
    dot_size_df = dot_size_df.stack().reset_index()
    dot_size_df.columns = ['Category', 'Gene', 'Size']


    dot_color_df = (df.groupby("group_splits", observed=True)[genes].mean()).stack().reset_index()
    dot_color_df.columns = ['Category', 'Gene', 'Color']

    merged_df = pd.merge(dot_size_df, dot_color_df, on=['Category', 'Gene'])


    if swap is True:
        x_cat = "Gene"
        y_cat = "Category"
    else:
        x_cat = "Category"
        y_cat = "Gene"
        
    
    fig = px.scatter(
        merged_df,
        x=x_cat,
        y=y_cat,
        size='Size',
        color='Color',
        color_continuous_scale = ["blue","white","red"],
        labels = {"Color":"Mean Expression In Group", "Size" : "Fraction of cells in group"},
        size_max=15
    )
    
    fig.update_layout(
        xaxis_title_text='',
        yaxis_title_text='',
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )

    
    return fig


# HEATMAP PLOT - SNS CLUSTERMAP CODE
def plot_heatmap(genes, category, category_subset, transform):

    category = sh.categoricals_dict[category]

    genes = list(set(genes))

    if transform == "TPM":
        df = pd.read_feather("./quipi_data/quipi_raw_tpm.feather", columns=[category] + list(genes))
    elif transform == "Log2(TPM)":
        df = pd.read_feather("./quipi_data/quipi_log2_tpm.feather", columns=[category] + list(genes))

    df = df[df[category].isin(category_subset)]
    df = df.groupby(category)[genes].median().T

    df = df[df.std(axis=1) > 0]


    fig = sns.clustermap(df, vmin=-2, vmax=2, cmap='coolwarm',
                         z_score=0,
                         xticklabels=True, yticklabels=True,
                         #figsize=(10, max(20, 1.25 * len(genes))),
                         dendrogram_ratio=(0.1, 0.05))

    fig.ax_heatmap.set_xlabel('')
    fig.ax_heatmap.set_ylabel('')
    fig.tick_params(axis='y', rotation=0)

    plt.subplots_adjust(
        left=.01,
        right=.9,
        top=1,
        bottom=.05
    )

    cbar_pos = [0.025, 0.95, 0.025, 0.025]  # [left, bottom, width, height]
    fig.ax_cbar.set_position(cbar_pos)
    fig.ax_cbar.set_title('Z-Score', fontdict={"fontsize":10})
    fig.ax_heatmap.set_yticklabels(fig.ax_heatmap.get_yticklabels(), size = 8)

    fig

    return fig

# RIDGELINE PLOT - CURRENTLY NOT IMPLEMENTED
def plot_ridgeline():

    sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

    df = pd.read_feather("/Users/hwismer/Documents/QuIPI/quipi_data/quipi_log2_tpm.feather")

    genes = ["C1QA", "FOXP3", "FOS", "JUNB", "CXCL3", "NCAM1", "DLK1", "NTRK3", "IGF2", "PAPPA2",]
             #'RAB33A','FFAR4','SERPINF2','GLYAT','SYT6','LINC02694','POU3F3','EMX1','LCP2','GBP3','APOBEC3H','CACYBP','CSF1',
             #'TNFSF14','CCL4L2','BPGM','C3','DNAJB4','CD200R1','ULBP1']
    #genes = ["C1QA", "FOXP3", "JUNB"]
    category = "archetype"
    sub_df = df[genes + [category]]

    tidy = pd.melt(sub_df,
          id_vars=category,
          value_vars=genes,
          var_name='Gene',
          value_name='Expression')
    

    g = sns.FacetGrid(tidy, row = "Gene", hue=category, aspect=1, height=3)

    g.map(sns.kdeplot, "Expression",
        bw_adjust=1, clip_on=False, 
        fill=True, alpha=.5, linewidth=1.5)

    #g.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)

    for ax, row_label in zip(g.axes.flat, g.row_names):
        ax.text(0, .2, row_label, fontweight="bold", color=None,
                ha="left", va="center", transform=ax.transAxes)
        
    g.figure.subplots_adjust(hspace=-.25)
    #g.figure.subplots_adjust(top=1.1, bottom=.05)

    g.add_legend()
    g.set_titles("")
    g.set(yticks=[], ylabel="")
    g.despine(bottom=False, left=True)
    #g.figure.set_size_inches(fixed_aspect * fixed_height, len(genes) * 5)
    #plt.tight_layout()
    #plt.subplots_adjust(right=.9)
    g.add_legend(bbox_to_anchor=(1, 0.5))

    return g