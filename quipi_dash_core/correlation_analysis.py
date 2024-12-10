import shared as sh
import numpy as np
import pandas as pd
import plotly.express as px


def gene_correlation_heatmap(genes, indications, method, compartments, archetypes, tissues, transform):

        if transform == "Raw":
            input_arr = pd.read_feather("./data/quipi_raw_tpm.feather", columns=sh.non_genes + list(genes))
        elif transform == "Log2":
            input_arr = pd.read_feather("./data/quipi_log2_tpm.feather", columns==sh.non_genes + list(genes))

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