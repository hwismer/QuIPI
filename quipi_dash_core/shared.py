from pathlib import Path

import pandas as pd

app_dir = Path(__file__).parent

quipi_raw = pd.read_pickle("./data/clean/quipi_raw_tpm.pi")
quipi_log10 = pd.read_pickle("./data/clean/quipi_log10_tpm.pi")
quipi_log2 = pd.read_pickle("./data/clean/quipi_log2_tpm.pi")
quipi_flow = pd.read_pickle("./data/clean/flow_mat.pi")

pancan_only_raw = quipi_raw[quipi_raw["archetype"] != "Unclassified"]

categoricals = ["patient", "sample_name",
             "indication", "sample_type","sample_type_cat",
             "compartment", "archetype"]

non_genes = ["patient", "sample_name",
             "indication", "sample_type","sample_type_cat",
             "compartment", "archetype", 
             "x_umap1", "x_umap2"]

genes = list(set(quipi_raw.columns) - set(non_genes))

#genes = ["DLK1", "NCAM1", "IGF2", "PLAG1","LY6H", "MDK", "NTRK3", "FGFR1", "NTRK3", "SLC7A3"]

categoricals_dict = {"Patient" : "patient",
                     "Indication" : "indication", 
                     "Tissue" : "sample_type_cat",
                     "Compartment" : "compartment", 
                     "Archetype" : "archetype"}

tissue_dict = {"Tumor" : "T",
               "Normal" : "N"}

transformations = {"Raw" : quipi_raw,
                   "Log2" : quipi_log2,
                   "Log10" : quipi_log10}

corr_methods = {"Spearman" : "spearman",
                "Pearson" : "pearson"}

flow_scores = {"Myeloid" : "Myelo_score",
               "T Cell" : 'T_score',
               "Stroma" : "Stroma_score",
               "T Reg" : "Treg_score",
               "CD4" : "CD4_score",
               "CD8" : "CD8_score",
               "Macrophage" : "Mac_score",
               "Monocyte" : "Mono_score",
               "cDC1" : "cDC1_score",
               "cDC2" : "cDC2_score",
               "Exhausted" : "Ex_score"}

colors_pancan = {
    'IR CD8 Mac' : '#ed1e21',
    'IR CD8 Mono' : '#f06ba8',
    'IR CD4 Mac' : '#7e1515',
    'IS CD8' : '#128042',
    'IS CD4' : '#98ca3a',
    'TC Mac' : '#2c276b',
    'TC DC' : '#4a87c7',
    'MC DC2' : '#7f7f7f',
    'MC DC1' : 'black',
    'ID CD4 Mac' : '#fdd80d',
    'ID Mono' : '#b8882c',
    'ID CD8 Mac' : '#f68c20'}

colors_indic = {
    'BLAD' : 'gold',
    'CRC' : 'darkturquoise',
    'GBM' : 'grey50',
    'GYN' : 'salmon1', 
    'HEP' : 'cornflowerblue',
    'HNSC' : 'olivedrab',
    'KID' : '#ed1e21',
    'LUNG' : '#2c276b',
    'MEL' : 'gold4',
    'PDAC' : 'yellow3',
    'PNET' : 'plum4', 
    'SRC' : 'hotpink'}