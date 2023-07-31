from shiny import *

import numpy as np
#If uncomment below, plot will not be opened in seperate window
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import os

app_ui = ui.page_fluid(
    ui.input_slider("n", "Slider", 0, 100, 20),
    ui.output_plot("myplot"),
)

current_directory = os.getcwd()
relative_path = os.path.join(current_directory, "myapp", "venv", "bin", "eg.xlsx")
df = pd.read_excel(relative_path)

num_rows = len(df)


def server(input: Inputs, output: Outputs, session: Session):
    @output
    @render.plot(alt="A bar")

    def myplot():
        NumAA = list(range(1, num_rows + 1))
        Values = df["average"].tolist()
        
        sd = df["SD"].tolist()
        LABELS = df[2].tolist()
        
        plt.title("Predicted HDX rates")
        plt.xlabel("Peptide Fragment")
        plt.ylabel("Mean of Models")

        plt.bar(NumAA, Values, align='center')

        plt.xticks(NumAA, LABELS)
        plt.xticks(rotation=30, ha='right')
        plt.yticks(np.arange(0, 1, step=0.025), minor=True)

        plt.errorbar(NumAA, Values, yerr = sd, fmt="o", color="k")

        plt.tight_layout()
        
        plt.show()

app = App(app_ui, server, debug=True)
