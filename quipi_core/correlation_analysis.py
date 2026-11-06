import quipi_shared as sh
import numpy as np
import pandas as pd
import plotly.express as px
import scipy
import gene_factor as gf


# Main correlation visualization function
# User specifies all the categories they want to look in. Simplest correlation function.
def gene_correlation_heatmap(genes, indications, method, compartments, archetypes, tissues, transform):

    if transform == "TPM":
        input_arr = pd.read_feather("./quipi_data/quipi_raw_tpm.feather", columns=sh.non_genes + list(genes))
    elif transform == "Log2(TPM)":
        input_arr = pd.read_feather("./quipi_data/quipi_log2_tpm.feather", columns=sh.non_genes + list(genes))

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

    fig = px.imshow(corr_lower_tri.astype(object).fillna(""),
                    zmin = -1, zmax =1,
                    color_continuous_scale = "RdBu_r",
                    text_auto=True)
    
    fig.update_layout(template='simple_white',autosize=True)
    fig.update_xaxes(showgrid=False, showline=False)
    fig.update_yaxes(showgrid=False, showline=False)

    return fig


#
def compartment_correlation_heatmap(genes, compartment1, compartment2, transform, method, indications, tissues, archetypes):

    if transform == "TPM":
        input_arr = pd.read_feather("./quipi_data/quipi_raw_tpm.feather", columns=sh.non_genes + list(genes))
    elif transform == "Log2(TPM)":
        input_arr = pd.read_feather("./quipi_data/quipi_log2_tpm.feather", columns=sh.non_genes + list(genes))

    input_arr = input_arr[input_arr["indication"].isin(indications)]
    input_arr = input_arr[input_arr["archetype"].isin(archetypes)]
    input_arr = input_arr[input_arr["sample_type_cat"].isin(tissues)]

    comp1 = input_arr[input_arr["compartment"] == compartment1]
    comp2 = input_arr[input_arr["compartment"] == compartment2]

    comp1 = comp1[np.logical_not(comp1.duplicated(subset="sample_name",keep="last"))]
    comp2 = comp2[np.logical_not(comp2.duplicated(subset="sample_name",keep="last"))]

    comp1_samples = comp1["sample_name"]
    comp2_samples = comp2["sample_name"]

    merged_samples = pd.merge(comp1_samples, comp2_samples, on=['sample_name'], how='inner')

    comp1_sub = comp1[comp1["sample_name"].isin(merged_samples["sample_name"])].sort_values("sample_name").reset_index(drop=True)
    comp2_sub = comp2[comp2["sample_name"].isin(merged_samples["sample_name"])].sort_values("sample_name").reset_index(drop=True)

    if (comp1_sub["sample_name"] == comp2_sub["sample_name"]).all():
        comp1_sub_genes = comp1_sub[genes]
        comp2_sub_genes = comp2_sub[genes]
        
        df = pd.DataFrame()
        feat1s = []
        feat2s = []
        corrs = []
        p_values = []
        
        for gene1 in genes:
            for gene2 in genes:
                feat1s.append(gene1)
                feat2s.append(gene2)
                if method == "Spearman":
                    corr, p_value = scipy.stats.spearmanr(comp1_sub_genes[gene1], comp2_sub_genes[gene2])
                elif method == "Pearson":
                    corr, p_value = scipy.stats.pearsonr(comp1_sub_genes[gene1], comp2_sub_genes[gene2])
                corrs.append(corr)
                p_values.append(p_value)
        
        df['Feature_1'] = feat1s
        df['Feature_2'] = feat2s
        df['Correlation'] = corrs
        df['p_value'] = p_values

        corr_df = df.pivot(index='Feature_1', columns='Feature_2', values="Correlation")
        mask = np.logical_not(np.tri(N=len(corr_df),dtype=bool))
        corr_df = corr_df.mask(mask)

        corr_df = corr_df.dropna(axis=0, how = "all")
        corr_df = corr_df.dropna(axis=1, how = "all")
        
        fig = px.imshow(corr_df.astype(object).fillna(""),
                zmin = -1, zmax =1,
                color_continuous_scale = "RdBu_r",
                text_auto=True,
                )
        
        fig.update_layout(template='simple_white',autosize=True,xaxis_title="",yaxis_title="")
        fig.update_xaxes(showgrid=False, showline=False)
        fig.update_yaxes(showgrid=False, showline=False)

        return fig


# Unused
# def gene_corr_df():
#     genes = list(input.corr_gene_input())
#     indications = input.corr_indication()
#     method = sh.corr_methods[input.corr_method_input()]
#     compartments = input.corr_compartment()
#     archetypes = input.corr_archetype()
#     tissues = [sh.tissue_dict[tis] for tis in input.corr_tissue()]

#     transform = input.corr_transform()

#     input_arr = sh.transformations[transform]
#     input_arr = input_arr[input_arr["indication"].isin(indications)]
#     input_arr = input_arr[input_arr["compartment"].isin(compartments)]
#     input_arr = input_arr[input_arr["archetype"].isin(archetypes)]
#     input_arr = input_arr[input_arr["sample_type_cat"].isin(tissues)]
#     input_arr = input_arr[genes]

#     return input_arr




##### Categorical correlation table (currently not implemented).
# User chooses a category (compartment, archetype, etc) and correlations calculated for a gene within that category.
# def categorical_correlation_table(gene,category,categories,range,progress):

#     log2_df = pd.read_feather("./quipi_data/quipi_log2_tpm.feather")
#     log2_df = log2_df[log2_df[sh.categoricals_dict[category]].isin(categories)][sh.genes]

#     df = pd.DataFrame()
#     genes,corrs,p_values = [],[],[]

#     for count, gene2 in enumerate(log2_df.columns):
#         if gene != gene2:
#             corr, p_value = scipy.stats.spearmanr(log2_df[gene].astype("float32"), log2_df[gene2].astype("float32"), nan_policy="omit")
#             if corr >= range[0] and corr <= range[1]:
#                 genes.append(gene2)
#                 corrs.append(corr)
#                 p_values.append(p_value)
#         progress.set(count, message = "Calculating")

#     df['Gene'] = genes
#     df['Spearman R'] = corrs
#     df['P-Value'] = p_values

            
#     return df.sort_values(["Spearman R"], ascending=False) 


# Calculates correlations for specified genes in compartment1 with all genes in each compartment of compartment2
# Only does compartment right now, but could potentially integrate with the categorical correlation function above?
def cross_compartment_correlation_table(genes, compartment1, compartment2, range, transform, method, progress):

    if transform == "TPM":
        input_arr = pd.read_feather("./quipi_data/quipi_raw_tpm.feather")
    elif transform == "Log2(TPM)":
        input_arr = pd.read_feather("./quipi_data/quipi_log2_tpm.feather")

    comp1 = input_arr[input_arr["compartment"] == compartment1]
    comp1 = comp1[np.logical_not(comp1.duplicated(subset="sample_name",keep="last"))]
    comp1_samples = comp1["sample_name"]

    genes1 = []
    compartments1 = []
    genes2 = []
    compartments2 = []
    corrs = []
    p_values = []


    count = 0
    for compartment in compartment2:
        
        comp2 = input_arr[input_arr["compartment"] == compartment]
        comp2 = comp2[np.logical_not(comp2.duplicated(subset="sample_name",keep="last"))]
        comp2_samples = comp2["sample_name"]
        
        merged_samples = pd.merge(comp1_samples, comp2_samples, on=['sample_name'], how='inner')
        
        comp1_sub = comp1[comp1["sample_name"].isin(merged_samples["sample_name"])].sort_values("sample_name").reset_index(drop=True)
        comp2_sub = comp2[comp2["sample_name"].isin(merged_samples["sample_name"])].sort_values("sample_name").reset_index(drop=True)
        
        if (comp1_sub["sample_name"] == comp2_sub["sample_name"]).all():
            for gene1 in genes:
                for gene2 in comp2_sub[sh.genes].columns:
                    if method == "Spearman":
                        corr, p_value = scipy.stats.spearmanr(comp1_sub[gene1], comp2_sub[gene2])
                    elif method == "Pearson":
                        corr, p_value = scipy.stats.pearsonr(comp1_sub[gene1], comp2_sub[gene2])
                    if corr >= range[0] and corr <= range[1]:
                        genes1.append(gene1)
                        compartments1.append(compartment1)
                        genes2.append(gene2)
                        compartments2.append(compartment)
                        corrs.append(corr)
                        p_values.append(p_value)
                    count += 1
                    progress.set(count, message="Calculating - I'm Accurate!", detail="This could take a while.")

    df = pd.DataFrame()
    df['Gene1'] = genes1
    df["Compartment1"] = compartments1
    df["Gene2"] = genes2
    df["Compartment2"] = compartments2
    df[method + ' R'] = corrs
    df['P-Value'] = p_values

    return df .sort_values([method + " R"], ascending=False)
        



# Correlates two values with the chosen method, returning (R, p-value)
def correlate(x, y, method):
    if method == "Spearman":
        return scipy.stats.spearmanr(x, y)
    elif method == "Pearson":
        return scipy.stats.pearsonr(x, y)


# Correlates a gene-signature score in compartment A with a gene-signature score in compartment B.
# Signature scores use all samples in each compartment: per-gene Z-Score within the compartment,
# averaged per sample.
# Returns a scatter plot of the two scores and a table of the signature and per-gene correlations.
def signature_cross_compartment_correlation(genes_a, compartment_a, genes_b, compartment_b, method, color_by):

    genes_a, genes_b = list(genes_a), list(genes_b)
    scores_a = gf.calculate_gene_factor_score_all_patients(genes_a, [compartment_a])
    scores_b = gf.calculate_gene_factor_score_all_patients(genes_b, [compartment_b])

    # Keep one row per sample in each compartment, then pair samples present in both
    scores_a = scores_a.drop_duplicates(subset="sample_name", keep="last").set_index("sample_name")
    scores_b = scores_b.drop_duplicates(subset="sample_name", keep="last").set_index("sample_name")
    shared_samples = scores_a.index.intersection(scores_b.index)
    scores_a = scores_a.loc[shared_samples]
    scores_b = scores_b.loc[shared_samples]

    # Table: the signature correlation first, then every gene in A against every gene in B
    pairs = [("Signature A", "factor_score", "Signature B", "factor_score")]
    pairs += [(gene_a, gene_a, gene_b, gene_b) for gene_a in genes_a for gene_b in genes_b]

    rows = []
    for name_a, col_a, name_b, col_b in pairs:
        r, p_value = correlate(scores_a[col_a], scores_b[col_b], method)
        rows.append([name_a, compartment_a, name_b, compartment_b, round(r, 3), p_value])

    table = pd.DataFrame(rows, columns=["Feature1", "Compartment1", "Feature2", "Compartment2",
                                        method + " R", "P-Value"])
    r, p_value = rows[0][4], rows[0][5]

    # Scatter plot of the two signature scores
    plot_df = scores_a[["indication", "archetype"]].copy()
    plot_df["score_a"] = scores_a["factor_score"]
    plot_df["score_b"] = scores_b["factor_score"]
    plot_df = plot_df.reset_index()

    color_col = sh.categoricals_dict[color_by]

    fig = px.scatter(plot_df, x="score_a", y="score_b",
                     color=color_col, color_discrete_map=sh.color_dict[color_col],
                     category_orders=sh.quipi_orders,
                     hover_data=["sample_name", "indication", "archetype"],
                     labels={"score_a": f"{compartment_a} Signature Score",
                             "score_b": f"{compartment_b} Signature Score",
                             color_col: color_by},
                     title=f"{method} R = {r:.3f}, p = {p_value:.2e}, n = {len(plot_df)}")

    fig.update_traces(marker=dict(size=10))
    fig.update_layout(template="simple_white", autosize=True)

    return fig, table
