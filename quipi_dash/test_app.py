import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from shiny import App, render, ui

# Load a sample dataset
df = sns.load_dataset("iris")

# Define the UI
app_ui = ui.page_fluid(
    ui.input_selectize("species", "Choose Species", choices=list(df['species'].unique()), multiple=True),
    ui.output_ui("plots")
)

# Define the server logic
def server(input, output, session):
    @output
    @render.ui
    def plots():
        selected_species = input.species()
        
        if not selected_species:
            return ui.div("Please select at least one species.")
        
        plot_outputs = []
        for species in selected_species:
            fig, ax = plt.subplots()
            sns.histplot(df[df['species'] == species]['sepal_length'], kde=True, ax=ax)
            ax.set_title(f"Sepal Length Distribution for {species}")
            plot_id = f"plot_{species}"
            plt.tight_layout()
            output_plot = ui.output_plot(plot_id)
            plot_outputs.append(ui.div(output_plot))
            
            # Dynamically create plots
            @output
            @render.plot(alt="Distribution plot")
            def _(plot_id=plot_id, species=species):
                return fig
        
        return ui.div(*plot_outputs)

# Create the Shiny app
app = App(app_ui, server)

# Run the app
#app.run()