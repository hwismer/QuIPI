import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib
import humu_shared as hsh
import plotly.graph_objects as go

def plot_sc_box(gene, x_cat, x_cat_subset, groupby, splitby, sample_aggr):
    cols = {gene, x_cat, groupby, splitby, "Mouse"} - {"---"}
    input_arr = pd.read_feather("./quipi_humu_data/quipi_humu_adata_clean_full_PROC.feather", columns = tuple(cols))
    input_arr = input_arr[input_arr[x_cat].isin(x_cat_subset)]

    if sample_aggr:
        cats = list(set([x_cat, groupby, splitby]))
        keep_cats = []
        for cat in cats:
            if cat != "---":
                 keep_cats.append(cat)
                 

        input_arr = input_arr.groupby(keep_cats + ["Mouse"])[gene].mean().reset_index()

    splitby = splitby if splitby != "---" else None
    groupby = groupby if groupby != "---" else x_cat

    colors = hsh.groupby_colors[groupby]
    orders = hsh.humu_orders[groupby]


    fig = px.box(input_arr, x = x_cat, y = gene, color = groupby, facet_col=splitby, 
                    facet_col_wrap=2,
                    points = "outliers",
                    facet_col_spacing=0.001,
                    facet_row_spacing=0.03,
                    color_discrete_map=colors,
                    category_orders={x_cat:orders}
                    )
    fig.update_traces(marker=dict(size=2), selector=dict(type='box'))

    return fig


def plot_dotplot(genes, groupby, groups, splitby, splits, swap):

    CUTOFF=0

    if groupby == splitby:
        splitby = "---"

    if splitby == "---":
        df = pd.read_feather("./quipi_humu_data/quipi_humu_adata_clean_full_PROC.feather", 
                         columns = genes + [groupby])
        df = df[df[groupby].isin(groups)]
        df["group_splits"] = df[groupby].astype(str)
        
    else:
        df = pd.read_feather("./quipi_humu_data/quipi_humu_adata_clean_full_PROC.feather", 
                         columns = genes + [groupby] + [splitby])
        df = df[(df[groupby].isin(groups)) & (df[splitby].isin(splits))]
        df["group_splits"] = df[groupby].astype(str) + '_' + df[splitby].astype(str)

    
    
    expressed = df[genes] > 0
    expressed["group_splits"] = df["group_splits"]

    dot_size_sum = expressed.groupby("group_splits")[genes].sum()
    dot_size_total = expressed.groupby("group_splits")[genes].count()
    
    dot_size_df = dot_size_sum / dot_size_total
    dot_size_df = dot_size_df.stack().reset_index()
    dot_size_df.columns = ['Category', 'Gene', 'Size']


    dot_color_df = (df.groupby("group_splits", observed=True)[genes].mean()).stack().reset_index()
    dot_color_df.columns = ['Category', 'Gene', 'Color']

    merged_df = pd.merge(dot_size_df, dot_color_df, on=['Category', 'Gene'])

    
    all_y_ticks = merged_df['Category'].unique()

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



def humu_box_comparison_human(human_gene, human_compartment_filters, human_transformation, human_groupby):

    if human_transformation == "TPM":
            input_arr = pd.read_feather("./quipi_humu_data/quipi_raw_tpm.feather", columns=hsh.quipi_cats + [human_gene])
    elif human_transformation == "Log2(TPM)":
            input_arr = pd.read_feather("./quipi_humu_data/quipi_log2_tpm.feather", columns=hsh.quipi_cats + [human_gene])
    
    input_arr = input_arr[input_arr["compartment"].isin(human_compartment_filters)]

    human_groupby = hsh.quipi_cats_dict[human_groupby]
    colors = hsh.quipi_colors[human_groupby]

    fig = px.box(input_arr, x = "compartment", y = human_gene, 
                 color = human_groupby, color_discrete_map=colors,
                 category_orders={"compartment" : hsh.compartment_order, "archetype" : hsh.archetype_order},
                 labels={"compartment" : "Compartment"})

    return fig


def humu_box_comparison_mouse(mouse_gene, mouse_compartment_filters, sample_aggr, mouse_groupby):

    input_arr = pd.read_feather("./quipi_humu_data/quipi_humu_adata_clean_full_PROC.feather", columns = [mouse_gene,"Mouse"] + hsh.categoricals_opts)
    input_arr = input_arr[input_arr["Compartment"].isin(mouse_compartment_filters)]
    groupby_cols = None

    if mouse_groupby == "---":
        mouse_groupby = "Compartment"


    if sample_aggr:
        if mouse_groupby == "Compartment":
             groupby_cols = ["Compartment", "Mouse"]
        else:
             groupby_cols = ["Mouse", "Compartment", mouse_groupby]

        input_arr = input_arr.groupby(groupby_cols)[mouse_gene].mean().reset_index()


    colors = hsh.groupby_colors[mouse_groupby]

    fig = px.box(input_arr, x = "Compartment", y = mouse_gene, 
                color = mouse_groupby, 
                color_discrete_map=colors,
                category_orders={"Compartment" : hsh.humu_compartments,})

    return fig

    

