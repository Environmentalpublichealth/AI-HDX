import pandas as pd
from time import time
import os

from shiny import App, Inputs, Outputs, Session, render, ui
from shiny.types import FileInfo

app_ui = ui.page_fluid(
    ui.row(  
        ui.column(
            6, #length of input box
            ui.layout_sidebar(
                ui.panel_sidebar(
                    ui.input_file("file1", "Choose CSV File", accept=[".csv"], multiple=False),
                    width=5
                ),
                ui.panel_main(
                    ui.output_ui("contents1"),
                ),
            ),
        ),
        ui.column(
            6, #length of input box
            ui.layout_sidebar(
                ui.panel_sidebar(
                    ui.input_file("file2", "Choose TXT File", accept=[".txt"], multiple=False),
                    width=5
                ),
                ui.panel_main(
                    ui.output_ui("contents2"),
                ),
            ),
            int = 3 #space between input boxes
        ),
    ),
)

def server(input: Inputs, output: Outputs, session: Session):
    @output
    @render.ui
    def contents1():
        if input.file1() is None:
            return "Please upload a csv file"
        f: list[FileInfo] = input.file1()
        df = pd.read_csv(f[0]["datapath"])
        
        # save the data to disk
        filename_csv = "zig.kuang" + str(time()) + ".csv"
        df.to_csv(os.path.join("C:\\Users\\zigku\\OneDrive\\Desktop\\Anna", filename_csv))
        
        #save file to variable
        df_csv = pd.read_csv(os.path.join("C:\\Users\\zigku\\OneDrive\\Desktop\\Anna", filename_csv))

        # return visible table with data
        return ui.HTML(df.to_html(classes="table table-striped"))

    @output
    @render.ui
    
    def contents2():
        if input.file2() is None:
            return "Please upload a txt file"
        f: list[FileInfo] = input.file2()
        df = pd.read_csv(f[0]["datapath"], delimiter="\t")
        
        # save the data to disk
        filename_txt = "zig.kuang" + str(time()) + ".txt"
        df.to_csv(os.path.join("C:\\Users\\zigku\\OneDrive\\Desktop\\Anna", filename_txt))
        
        #save file to variable
        df_txt = pd.read_csv(os.path.join("C:\\Users\\zigku\\OneDrive\\Desktop\\Anna", filename_txt), delimiter="\t")

        # return visible table with data
        return ui.HTML(df.to_html(classes="table table-striped"))

app = App(app_ui, server)


'''
app_ui = ui.page_fluid(
    ui.row(  
        
        ui.column(
            6,
            ui.layout_sidebar(
                ui.panel_sidebar(
                    ui.input_file("file1", "Choose CSV File", accept=[".csv"], multiple=False),
                    ui.input_checkbox("header", "Header", True),
                    width=5
                ),
                ui.panel_main(
                    ui.output_ui("contents1"),

                ),
            ),
            
        ),
        ui.column(
            6,
            ui.layout_sidebar(
                ui.panel_sidebar(
                    ui.input_file("file2", "Choose TXT File", accept=[".txt"], multiple=False),
                    ui.input_checkbox("header", "Header", True),
                    width=5
                ),
                ui.panel_main(
                    ui.output_ui("contents2"),

                ),
            ),
            int = 3
        ),
    ),
)


def server(input: Inputs, output: Outputs, session: Session):
    @output
    @render.ui
    def contents1():
        if input.file1() is None:
            return "Please upload a csv file"
        f: list[FileInfo] = input.file1()
        df = pd.read_csv(f[0]["datapath"], header=0 if input.header() else None)
        # save the data to disk
        # os.makedirs('AI-HDX_app/save_file', exist_ok=True)  
        # df.to_csv('AI-HDX_app/save_file/upload.csv') 
        return ui.HTML(df.to_html(classes="table table-striped"))
    def contents2():
        if input.file2() is None:
            return "Please upload a txt file"
        g: list[FileInfo] = input.file2()
        df = pd.read_csv(g[0]["datapath"], header=0 if input.header() else None)
        # save the data to disk
        # os.makedirs('AI-HDX_app/save_file', exist_ok=True)  
        # df.to_csv('AI-HDX_app/save_file/upload.csv') 
        return ui.HTML(df.to_html(classes="table table-striped"))


app = App(app_ui, server)
'''