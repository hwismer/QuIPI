import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import csv


# All categorial columns for the underlying data minus the UMAP coordinates

# Everything that isn't a gene, including umap coordinates.
non_genes = ["patient", "sample_name",
             "indication", "sample_type","sample_type_cat",
             "compartment", "archetype", 
             "x_umap1", "x_umap2"]

categoricals = ["patient", "sample_name", "indication", "sample_type","sample_type_cat", "compartment", "archetype"]
indications = ['LUNG', 'HEP', 'ADR', 'GBM', 'CRC', 'BRC', 'KID', 'MEL', 'PNET', 'GYN', 'HNSC', 'SI', 'SRC', 'GALL', 'PDAC', 'BLAD']
compartments = ['Live', 'Tumor', 'T-Reg', 'Myeloid', 'Stroma', 'T Cell']
archetypes = ['Unclassified', 'IR CD8 Mac', 'IR CD8 Mono', 'IR CD4 Mac', 'IS CD8','IS CD4', 'TC Mac','TC DC','MC DC2','MC DC1','ID CD4 Mac','ID Mono','ID CD8 Mac']

categorical_choices = {"Compartment":compartments, 
                       "Archetype": archetypes, 
                       "Indication": indications}


##### ORDERS
compartment_order = ["Tumor", "T Cell", "T-Reg", "Myeloid", "Stroma", "Live"]
indication_order = ['LUNG', 'HEP', 'ADR', 'GBM', 'CRC', 'BRC', 'KID', 'MEL', 'PNET', 'GYN', 'HNSC', 'SI', 'SRC', 'GALL', 'PDAC', 'BLAD']
archetype_order = ['Unclassified', 'IR CD8 Mac', 'IR CD8 Mono', 'IR CD4 Mac', 'IS CD8','IS CD4', 'TC Mac','TC DC','MC DC2','MC DC1','ID CD4 Mac','ID Mono','ID CD8 Mac']
sample_type_order = ["T1", "T2", "N1", "N2"]

quipi_orders = {"indication" : indication_order,
               "sample_type" : sample_type_order,
               "archetype" : archetype_order,
               "compartment" : compartment_order}

##### COLORS

indic_to_color = {'LUNG':'rgb(15, 8, 58)',
                  'HEP':'rgb(55, 77, 161)', 
                  'ADR':'rgb(248, 156, 116)', ####
                  'GBM':'rgb(69, 70, 75)', 
                  'CRC':'rgb(74, 135, 145)',
                  'BRC':'rgb(158, 185, 243)', ####
                  'KID':'rgb(124, 20, 20)',
                  'MEL':'rgb(15, 2, 2)',
                  'PNET':'rgb(82, 26, 113)',
                  'GYN':'rgb(153, 99, 87)',
                  'HNSC':'rgb(39, 111, 36)',
                  'SI':'rgb(220, 176, 242)', ####
                  'SRC':'rgb(281, 23, 110)',
                  'GALL':'rgb(246, 207, 113)', ####
                  'PDAC':'rgb(165, 171, 44)',
                  'BLAD':'rgb(212, 160, 64)'}

compartment_colors = {"T Cell": "#009292",
                      "Myeloid" : "#FF6DB6",
                      "Stroma" : "#006DDB",
                      "Tumor" : "#DB6D00",
                      "T-Reg" : "#B66DFF",
                      "Live" : "#004949"}



# Open the file in read mode
with open("./quipi_data/quipi_cols.csv", 'r') as file:
    reader = csv.reader(file)
    # Read the first row only (header)
    quipi_all_columns= next(reader)

with open("./quipi_data/quipi_flow_score_cols.csv", 'r') as file:
    reader = csv.reader(file)
    # Read the first row only (header)
    quipi_flow_columns= next(reader)

cats = ('patient','sample_name','indication', 'archetype', 'sample_type', 'sample_type_cat', 'x_umap1', 'x_umap2')
non_cats = tuple(set(quipi_flow_columns) - set(cats))
genes = list(set(quipi_all_columns) - set(non_genes))

categorical_data = pd.read_feather("./quipi_data/quipi_raw_tpm.feather", columns=cats)

# Cleaned categories for user mapped to underlying data column names
categoricals_dict = {"Indication" : "indication", 
                     "Sample Type" : "sample_type",
                     "Compartment" : "compartment", 
                     "Archetype" : "archetype"}

categorials_opts_dict = {"Indication" : indications,
                         "Sample Type" : ["T1", "T2", "N1", "N2"],
                         "Compartment" : compartments,
                         "Archetype" : archetypes}

categoricals_dict_reversed = {y:x for x,y in categoricals_dict.items()}

# Mapped names for each simple T/N indication
tissue_dict = {"Tumor" : "T",
               "Normal" : "N"}

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
    'Unclassified' : '#d6d893',
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
    'ADR' : ["Adrenal"], 
    'GBM' : ["Glioblastoma"], 
    'CRC' : ["Colorectal"], 
    'BRC' : ["Breast"], 
    'KID' : ["Kidney"], 
    'MEL' : ["Melanoma"], 
    'PNET' : ["Primitive Neuro-Ectodermal"], 
    'GYN' : ["Gynecological"], 
    'HNSC' : ["Head & Neck Squamous Cell Carcinoma"], 
    'SI' : ["Small Intestinal"], 
    'SRC' : ["Sarcoma"], 
    'GALL' : ["Gallbladder"], 
    'PDAC' : ["Pancreatic Ductal Adenocarcinoma"], 
    'BLAD' : ["Bladder"]
}

cancer_glossary_df = pd.DataFrame.from_dict(cancer_glossary,
                                            orient = "index",
                                            columns = ["Elaborated"],)
cancer_glossary_df["Abbreviation"] = cancer_glossary_df.index

color_dict = {"indication" : indic_to_color,
              "sample_type" : None,
              "compartment" : compartment_colors,
              "archetype" : colors_pancan
              }


def plot_cancer_glossary_table():
    df = cancer_glossary_df.sort_values("Abbreviation")
    colors = [indic_to_color[ind] for ind in df["Abbreviation"]]

    fig = go.Figure(data=[
        go.Table(
            cells=dict(values=[df[col] for col in df.columns],
                    fill_color=[colors],  # Apply row colors
                    align='center',
                    height=28,
                    line=dict(color='white', width=2),
                    font = dict(color = 'white', size = 16)))
    ])
    
    fig.update_layout(autosize=True,)
    return fig

def plot_indication_breakdown():
    ind_counts = categorical_data.groupby("patient")["indication"].unique().value_counts().reset_index()
    ind_counts["indication"] = [ind[0] for ind in ind_counts["indication"]]
    
    fig = px.bar(ind_counts, x = "indication", y = "count", 
                 color = "indication", color_discrete_map=indic_to_color)
    
    fig.update_layout(showlegend=False,xaxis_title="",yaxis_title="Count")
    fig.update_layout(title_x= .5, title_y = .98,font=dict(size=16))

    return fig

def plot_archetype_beakdown():
    arch_counts = categorical_data[["sample_name", "archetype"]].groupby("sample_name")["archetype"].unique().value_counts().reset_index()
    arch_counts["archetype"] = [ind[0] for ind in arch_counts["archetype"]]

    fig = px.bar(arch_counts, x = "archetype", y = "count", 
                 color = "archetype", color_discrete_map=colors_pancan)
    
    fig.update_layout(showlegend=False,xaxis_title="",yaxis_title="Count")
    fig.update_layout(title_x= .5, title_y = .98,font=dict(size=16))

    return fig
