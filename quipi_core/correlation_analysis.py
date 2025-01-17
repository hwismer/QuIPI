import shared as sh
import numpy as np
import pandas as pd
import plotly.express as px
import scipy

from shiny import ui
import asyncio


def gene_correlation_heatmap(genes, indications, method, compartments, archetypes, tissues, transform):

        if transform == "Raw":
            input_arr = pd.read_feather("./data/quipi_raw_tpm.feather", columns=sh.non_genes + list(genes))
        elif transform == "Log2":
            input_arr = pd.read_feather("./data/quipi_log2_tpm.feather", columns=sh.non_genes + list(genes))

        input_arr = input_arr[input_arr["indication"].isin(indications)]
        input_arr = input_arr[input_arr["compartment"].isin(compartments)]
        input_arr = input_arr[input_arr["archetype"].isin(archetypes)]
        input_arr = input_arr[input_arr["sample_type_cat"].isin(tissues)]
        input_arr = input_arr[list(genes)]

        corr_df = input_arr.corr(method=method)

        # Create a mask for the upper triangle
        mask = np.triu(np.ones_like(corr_df, dtype=bool))

        # Apply the mask to the correlation matrix
        corr_lower_tri = corr_df.mask(mask)

        corr_lower_tri = corr_lower_tri.dropna(axis=0, how = "all")
        corr_lower_tri = corr_lower_tri.dropna(axis=1, how = "all")

        fig = px.imshow(corr_lower_tri.fillna(""),
                        zmin = -1, zmax =1,
                        color_continuous_scale = "RdBu_r",
                        text_auto=True)
        
        fig.update_layout(template='simple_white',autosize=True)
        fig.update_xaxes(showgrid=False, showline=False)
        fig.update_yaxes(showgrid=False, showline=False)

        return fig


def gene_corr_df():
    genes = list(input.corr_gene_input())
    indications = input.corr_indication()
    method = sh.corr_methods[input.corr_method_input()]
    compartments = input.corr_compartment()
    archetypes = input.corr_archetype()
    tissues = [sh.tissue_dict[tis] for tis in input.corr_tissue()]

    transform = input.corr_transform()

    input_arr = sh.transformations[transform]
    input_arr = input_arr[input_arr["indication"].isin(indications)]
    input_arr = input_arr[input_arr["compartment"].isin(compartments)]
    input_arr = input_arr[input_arr["archetype"].isin(archetypes)]
    input_arr = input_arr[input_arr["sample_type_cat"].isin(tissues)]
    input_arr = input_arr[genes]

    return input_arr

def categorical_correlation_table(gene,category,categories,range,progress):


    log2_df = pd.read_feather("./data/quipi_log2_tpm.feather")
    log2_df = log2_df[log2_df[sh.categoricals_dict[category]].isin(categories)][sh.genes]

    df = pd.DataFrame()
    genes,corrs,p_values = [],[],[]

    for count, gene2 in enumerate(log2_df.columns):
        if gene != gene2:
            corr, p_value = scipy.stats.spearmanr(log2_df[gene].astype("float32"), log2_df[gene2].astype("float32"), nan_policy="omit")
            if corr >= range[0] and corr <= range[1]:
                genes.append(gene2)
                corrs.append(corr)
                p_values.append(p_value)
        progress.set(count, message = "Calculating")

    df['Gene'] = genes
    df['Spearman R'] = corrs
    df['P-Value'] = p_values

            
    return df.sort_values(["Spearman R"], ascending=False) 





    
    log2_df = pd.read_feather("./data/quipi_log2_tpm.feather")
    log2_df = log2_df[log2_df[sh.categoricals_dict[category]].isin(categories)][sh.genes]

    df = pd.DataFrame()
    genes,corrs,p_values = [],[],[]

    with ui.Progress(min=1, max = len(log2_df.columns)) as p:
        p.set(message="Calculating", detail="Please Wait")
        for count, gene2 in enumerate(log2_df.columns):
            if gene != gene2:
                corr, p_value = scipy.stats.spearmanr(log2_df[gene].astype("float32"), log2_df[gene2].astype("float32"), nan_policy="omit")
                if corr >= range[0] and corr <= range[1]:
                    genes.append(gene2)
                    corrs.append(corr)
                    p_values.append(p_value)

        if count % 15 == 0:
            p.set(count, message = "Processing")

    df['Gene'] = genes
    df['Spearman R'] = corrs
    df['P-Value'] = p_values

    
    return df.sort_values(["Spearman R"], ascending=False)
        
