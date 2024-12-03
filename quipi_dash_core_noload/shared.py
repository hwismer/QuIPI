from pathlib import Path
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import csv

#app_dir = Path(__file__).parent

#quipi_raw = pd.read_pickle("./data/clean/quipi_raw_tpm.pi")
#quipi_log10 = pd.read_pickle("./data/clean/quipi_log10_tpm.pi")
#quipi_log2 = pd.read_pickle("./data/clean/quipi_log2_tpm.pi")
#quipi_flow = pd.read_pickle("./data/clean/flow_mat.pi")



# TPM MATRICES

# All categorial columns for the underlying data minus the UMAP coordinates
categoricals = ("patient", "sample_name",
             "indication", "sample_type","sample_type_cat",
             "compartment", "archetype")

# Everything that isn't a gene, including umap coordinates.
non_genes = ("patient", "sample_name",
             "indication", "sample_type","sample_type_cat",
             "compartment", "archetype", 
             "x_umap1", "x_umap2")

indications = ('LUNG', 'HEP', 'ADR', 'GBM', 'CRC', 'BRC', 'KID', 'MEL', 'PNET', 'GYN', 'HNSC', 'SI', 'SRC', 'GALL', 'PDAC', 'BLAD')
compartments = ('Live', 'Tumor', 'Treg', 'Myeloid', 'Stroma', 'T_cell')
archetypes = ('Unclassified', 'ID Mono', 'ID CD8 Mac', 'IR CD8 Mono', 
              'MC DC2', 'IS CD8', 'TC Mac', 'TC DC', 'IR CD4 Mac', 
              'IR CD8 Mac', 'ID CD4 Mac', 'MC DC1', 'IS CD4')

with open("./quipi_raw_tpm.csv", 'r') as file:
    reader = csv.DictReader(file)
    quipi_all_columns = reader.fieldnames


# FLOW MATRIX

with open("./quipi_flow_scores.csv", 'r') as file:
    reader = csv.DictReader(file)
    quipi_flow_columns = reader.fieldnames

cats = ('patient','sample_name','indication', 'archetype', 'sample_type', 'sample_type_cat', 'x_umap1', 'x_umap2')
non_cats = tuple(set(quipi_flow_columns) - set(cats))


#genes = list(set(quipi_all_columns) - set(non_genes))



pancan_only_raw = pd.read_csv("./quipi_raw_tpm.csv", usecols=["archetype", "x_umap1", "x_umap2"])


genes = ["DLK1", "NCAM1", "IGF2", "PLAG1","LY6H", "MDK", "NTRK3", "FGFR1", "NTRK3", "SLC7A3"]


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
#transformations = {"Raw" : quipi_raw,
#                   "Log2" : quipi_log2,
#                   "Log10" : quipi_log10}

corr_methods = {"Spearman" : "spearman",
                "Pearson" : "pearson"}

# Mapped names for each feature score to its corresponding column name in the data
feature_scores = {"Myeloid Score" : "Myelo_score",
               "T Cell Score" : 'T_score',
               "Stroma Score" : "Stroma_score",
               "T Reg Score" : "Treg_score",
               "CD4 Score" : "CD4_score",
               "CD8 Score" : "CD8_score",
               "Macrophage Score" : "Mac_score",
               "Monocyte Score" : "Mono_score",
               "cDC1 Score" : "cDC1_score",
               "cDC2 Score" : "cDC2_score",
               "Exhaustion Score" : "Ex_score"}

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

cancer_glossary = {
    'LUNG' : ["Lung"], 
    'HEP' : ["Hepatobiliary"], 
    'ADR' : ["ADR"], 
    'GBM' : ["Glioblastoma"], 
    'CRC' : ["Colorectal"], 
    'BRC' : ["Breast"], 
    'KID' : ["Kidney"], 
    'MEL' : ["Melanoma"], 
    'PNET' : ["Primitive Neuro-Ectodermal"], 
    'GYN' : ["Gynecological"], 
    'HNSC' : ["Head & Neck Squamous Cell Carcinoma"], 
    'SI' : ["SI"], 
    'SRC' : ["Sarcoma"], 
    'GALL' : ["Gallbladder"], 
    'PDAC' : ["Pancreatic Ductal Adenocarcinoma"], 
    'BLAD' : ["Bladder"]
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
                    fill_color='pink',
                    font = dict(color = "black",size = 18),
                    align='center'),
        cells=dict(values=[df[col] for col in df.columns],
                #fill_color=[[colors_indic[color] for color in df["Abbreviation"]]],  # Apply row colors
                align='center',
                height=30,
                font = dict(color = 'black', size = 18)))
    ])
    
    fig.update_layout(autosize=False, width=600,height=400)
    return fig


def plot_indication_breakdown():
    ind_counts = pd.read_csv("./quipi_raw_tpm.csv", usecols=non_genes).groupby("patient")["indication"].unique().value_counts().reset_index()
    
    #quipi_raw.groupby("patient")["indication"].unique().value_counts().reset_index()
    ind_counts["indication"] = [ind[0] for ind in ind_counts["indication"]]
    #ind_counts["color"] = [colors_indic[indic] if indic in colors_indic else None for indic in ind_counts["indication"]]
    fig = px.bar(ind_counts, x = "indication", y = "count", color = "indication",color_discrete_sequence=px.colors.qualitative.Set1)
    
    
    fig.update_layout(showlegend=False,xaxis_title="",yaxis_title="n")
    fig.update_layout(autosize=False, width=800, height=450)
    fig.update_layout(title="n Cancer Indication", title_x= .5, title_y = .98,font=dict(size=16))

    return fig



def plot_archetype_beakdown():
    pass

def plot_():
    pass

print('here')