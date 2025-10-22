import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib
import humu_shared as hsh
import plotly.graph_objects as go
from scipy.stats import zscore
import seaborn as sns
from pandas.api.types import is_numeric_dtype

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


def plot_dotplot(genes, groupby, groups, splitby, splits, swap, scale):

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

    #merged_df["Scaled_Color"] = merged_df.groupby("Gene")["Color"].transform(zscore)


    
    all_y_ticks = merged_df['Category'].unique()

    if swap is True:
        x_cat = "Gene"
        y_cat = "Category"
    else:
        x_cat = "Category"
        y_cat = "Gene"
    
    if scale is True:
        merged_df["Scaled_Color"] = merged_df.groupby("Gene")["Color"].transform(zscore)
        color = "Scaled_Color"
        range_color = [-2.5,2.5]
    else:
        color = "Color"
        range_color=None
    
    fig = px.scatter(
        merged_df,
        x=x_cat,
        y=y_cat,
        size='Size',
        color=color,
        color_continuous_scale = ["blue","lightgray","red"],
        labels = {"Color":"Mean Expression In Group", 
                  "Size" : "Fraction of cells in group",
                  "Scaled_Color":"Z-Score"},
        size_max=15,
        range_color=range_color
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


# Plots a pancan umap for each gene in genes list within a certain compartment
def plot_humu_umap(genes, categories):
    
    
    input_arr = pd.read_feather("./quipi_humu_data/quipi_humu_adata_clean_full_PROC.feather", 
                                columns =  genes + categories + ["X_UMAP", "Y_UMAP"])
    

    combined_length = len(genes + categories)

    # Create main figure, size set by total number of plots needed
    n_col = 4
    n_rows = (combined_length + n_col - 1) // n_col
    horizontal_spacing = 0.05

    fig = make_subplots(rows =  n_rows , cols = n_col, 
                        subplot_titles=genes + categories,
                        vertical_spacing=.05,horizontal_spacing=horizontal_spacing,
                        shared_xaxes=True,shared_yaxes=True)

    # count will keep track of how many plots have been plotted to allow for calculation of correct row and column
    count = 0

    # Begin with the genes which should all have numeric values
    for gene in genes:
        
        row = count // n_col
        col = count % n_col

        # Calculations for correct colorbar placement (still a bit wonky)
        x_col_end = (col + 1) / n_col
        if col < n_col - 1:
            x_domain_end = x_col_end - (horizontal_spacing / 2)
        else:
            x_domain_end = x_col_end
        colorbar_x = x_domain_end - 0.002
        y_center = 1 - ((row + 0.5) / n_rows)


        scatter = go.Scattergl(
            x = input_arr["X_UMAP"], 
            y = input_arr["Y_UMAP"],
            mode = 'markers',
            marker=dict(
                size=4,
                color=input_arr[gene],
                colorscale='Viridis',
                # May need tweaking
                colorbar=dict(
                    x=colorbar_x,
                    y=y_center,
                    yanchor="middle",
                    lenmode="fraction",
                    len=0.75 / n_rows,
                    thickness=15, 
                    thicknessmode="pixels"
                ),
            )
        )

    
        fig.add_trace(scatter, row= row+1, col= col+1)
        count += 1

    # Proceed to non-gene options such as nFeature_RNA or Coarse Annotation
    # Mix of dtypes here needs consideration
    for cat in categories:

        row = count // n_col
        col = count % n_col

        # Check if the category is numeric, in which case give it a continuous color palette and a colorbar
        if is_numeric_dtype(input_arr[cat]):

            # Colorbar placement
            x_col_end = (col + 1) / n_col
            if col < n_col - 1:
                x_domain_end = x_col_end - (horizontal_spacing / 2)
            else:
                x_domain_end = x_col_end
            colorbar_x = x_domain_end - 0.002
            y_center = 1 - ((row + 0.5) / n_rows)


            scatter = go.Scattergl(x = input_arr["X_UMAP"], y = input_arr["Y_UMAP"],
                                   mode = 'markers',
                                   marker=dict(
                                       size=4,  # Adjust marker size if needed
                                       color=input_arr[cat],  # Color by the numerical value
                                       colorscale='Viridis',
                                       colorbar=dict(
                                            x=colorbar_x,
                                            y=y_center,
                                            yanchor="middle",
                                            lenmode="fraction", # Use 'len' as a fraction of total figure height
                                            len=0.75 / n_rows, # Set length relative to the subplot height (e.g., 75%)
                                            thickness=15, 
                                            thicknessmode="pixels"
                                        ),
                                    )
                                )
        # Otherwise if the category is qualitative, give it a discrete color palette
        else:

            color_sequence = px.colors.qualitative.Safe 
            unique_categories = input_arr[cat].unique()

            for i, category_value in enumerate(unique_categories):
                color = color_sequence[i % len(color_sequence)] 
                df_cat = input_arr[input_arr[cat].astype(str) == category_value]
                
                # Create the trace
                scatter = go.Scattergl(
                    x=df_cat["X_UMAP"], 
                    y=df_cat["Y_UMAP"],
                    mode='markers',
                    name=category_value, # Name for the legend
                    legendgroup=cat, 
                    text=category_value,
                    hoverinfo='text',
                    showlegend=True, 
                    marker=dict(
                        size=4,
                        color=color, # Assign single, discrete color
                    ),
                )
                fig.add_trace(scatter, row=row + 1, col=col + 1)

    
        fig.add_trace(scatter, row= row+1, col= col+1)
        count += 1
    
         

    fig.update_xaxes(scaleanchor="y", scaleratio=1, showticklabels=False)
    fig.update_yaxes(scaleanchor="x", scaleratio=1, showticklabels=False)
    fig.update_layout(height=300*n_rows, width=290*n_col, showlegend=False,
                     uirevision=True,)

    return fig





    

