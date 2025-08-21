import plotly.express as px
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np

import quipi_shared as sh
import gene_factor as gf

import seaborn as sns


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

def plot_heatmap(genes, category, category_subset, transform):

    category = sh.categoricals_dict[category]

    genes = list(set(genes))

    if transform == "TPM":
        df = pd.read_feather("./quipi_data/quipi_raw_tpm.feather", columns=[category] + list(genes))
    elif transform == "Log2(TPM)":
        df = pd.read_feather("./quipi_data/quipi_log2_tpm.feather", columns=[category] + list(genes))

    df = df[df[category].isin(category_subset)]

    fig = sns.clustermap(df.groupby(category)[genes].sum().T,
              z_score=0, vmin=-2, vmax=2,cmap='coolwarm')
    

    fig.ax_heatmap.set_yticks(np.arange(len(genes)) + 0.5)
    fig.ax_heatmap.set_yticklabels(genes)  

    fig.ax_heatmap.set_xlabel('')
    fig.ax_heatmap.set_ylabel('')
    fig.tick_params(axis='y', rotation=0)

    plt.subplots_adjust(
        left=.01,
        right=.9,
        top=.95,
        bottom=.2
    )

    cbar_pos = [0.01, 0.8, 0.05, 0.1]  # [left, bottom, width, height]
    fig.ax_cbar.set_position(cbar_pos)
    fig.ax_cbar.set_title('Z-Score', fontdict={"fontsize":10})

    #plt.tight_layout()

    return fig