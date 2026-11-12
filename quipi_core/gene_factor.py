import numpy as np
import pandas as pd
from scipy.stats import zscore
import quipi_shared as sh


# Z-scores each gene across the given samples and averages the Z-scores per sample
# into a "factor_score" column. Returns the categorical columns, gene Z-scores and score.
def add_signature_score(samples, gene_set):
    z_subset = samples[gene_set].apply(zscore)
    z_subset["factor_score"] = z_subset.mean(axis=1)
    return samples[sh.non_genes].merge(z_subset, left_index=True, right_index=True)


# Helper function to calculate a gene factor score in a compartment given a set of genes
def calculate_gene_factor_score(gene_set,compartment):

    input_arr = pd.read_feather("./quipi_data/quipi_log2_tpm.feather", columns=sh.non_genes + gene_set)
    keep = (input_arr["archetype"] != "Unclassified") & (input_arr["compartment"] == compartment)

    return add_signature_score(input_arr[keep], gene_set)


def calculate_gene_factor_score_all_patients(gene_set,compartments):

    input_arr = pd.read_feather("./quipi_data/quipi_log2_tpm.feather", columns=sh.non_genes + gene_set)
    input_arr = input_arr[input_arr["compartment"].isin(compartments)]

    return add_signature_score(input_arr, gene_set)


# Signature score within one compartment, using only samples from the chosen
# indications, archetypes and tissues (tissues as sample_type_cat values, e.g. "T", "N")
def calculate_gene_factor_score_subset(gene_set, compartment, indications, archetypes, tissues):

    input_arr = pd.read_feather("./quipi_data/quipi_log2_tpm.feather", columns=sh.non_genes + gene_set)
    keep = ((input_arr["compartment"] == compartment)
            & input_arr["indication"].isin(indications)
            & input_arr["archetype"].isin(archetypes)
            & input_arr["sample_type_cat"].isin(tissues))

    if not keep.any():
        raise ValueError(f"No {compartment} samples match the selected tissues, archetypes and indications.")

    return add_signature_score(input_arr[keep], gene_set)
