from pathlib import Path
import pandas as pd
import plotly.graph_objects as go

app_dir = Path(__file__).parent

quipi_raw = pd.read_pickle("./data/clean/quipi_raw_tpm.pi")
quipi_log10 = pd.read_pickle("./data/clean/quipi_log10_tpm.pi")
quipi_log2 = pd.read_pickle("./data/clean/quipi_log2_tpm.pi")
quipi_flow = pd.read_pickle("./data/clean/flow_mat.pi")

pancan_only_raw = quipi_raw[quipi_raw["archetype"] != "Unclassified"]

# All categorial columns for the underlying data minus the UMAP coordinates
categoricals = ["patient", "sample_name",
             "indication", "sample_type","sample_type_cat",
             "compartment", "archetype"]

# Everything that isn't a gene, including umap coordinates.
non_genes = ["patient", "sample_name",
             "indication", "sample_type","sample_type_cat",
             "compartment", "archetype", 
             "x_umap1", "x_umap2"]


genes = list(set(quipi_raw.columns) - set(non_genes))

#genes = ["DLK1", "NCAM1", "IGF2", "PLAG1","LY6H", "MDK", "NTRK3", "FGFR1", "NTRK3", "SLC7A3"]

# Cleaned categories for user mapped to underlying data column names
categoricals_dict = {"Patient" : "patient",
                     "Indication" : "indication", 
                     "Tissue" : "sample_type_cat",
                     "Compartment" : "compartment", 
                     "Archetype" : "archetype"}

categoricals_dict_reversed = {y:x for x,y in categoricals_dict.items()}

# Mapped names for each simple T/N indication
tissue_dict = {"Tumor" : "T",
               "Normal" : "N"}

# Mapped names for each data representation
transformations = {"Raw" : quipi_raw,
                   "Log2" : quipi_log2,
                   "Log10" : quipi_log10}

corr_methods = {"Spearman" : "spearman",
                "Pearson" : "pearson"}

# Mapped names for each feature score to its corresponding column name in the data
feature_scores = {"Myeloid" : "Myelo_score",
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

# Colors for each PanCan archetype
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

# Colors for each cancer indication
colors_indic = {
    'BLAD' : "#FFD700",
    'CRC' : "#00CED1",
    'GBM' : ' #7f7f7f',
    'GYN' : "#FF8C69", 
    'HEP' : "#6495ED",
    'HNSC' : "#6B8E23",
    'KID' : '#ed1e21',
    'LUNG' : '#2c276b',
    'MEL' : "#8B7500",
    'PDAC' : "#CDCD00",
    'PNET' : "#8B668B", 
    'SRC' : "#FF69B4"}

# Used to create abbreviation table on the front page.
cancer_glossary = {
    "BLAD" : ["Bladder"],
    "CRC" : ["Colorectal"],
    "GBM" : ["Glioblastoma"],
    "HEP" : ["Hepatobiliary"],
    "HNSC" : ["Head & Neck Squamous Cell Carcinoma"],
    "KID" : ["Kidney"],
    "LUNG" : ["Lung"],
    "MEL" : ["Melanoma"],
    "PDAC" : ["Pancreate Ductal Adenocarcinoma"],
    "PNET" : ["Primitive Neuro-Ectodermal"],
    "SRC" : ["Sarcoma"]
}
cancer_glossary_df = pd.DataFrame.from_dict(cancer_glossary,
                                            orient = "index",
                                            columns = ["Elaborated"],)
cancer_glossary_df["Abbreviation"] = cancer_glossary_df.index
cancer_glossary_df = cancer_glossary_df[["Abbreviation", "Elaborated"]]

def plot_cancer_glossary_table():
    df = cancer_glossary_df

    fig = go.Figure(data=[go.Table(
        header=dict(values=list(df.columns),
                    fill_color='white',
                    font = dict(color = "black",size = 18),
                    align='center'),
        cells=dict(values=[df[col] for col in df.columns],
                fill_color=[[colors_indic[color] for color in df["Abbreviation"]]],  # Apply row colors
                align='center',
                height=30,
                font = dict(color = 'white', size = 18)))
    ])
    
    fig.update_layout(autosize=False, width=600,height=800)
    return fig