import pandas as pd
import numpy as np
from time import time
import os
import glob
# import asyncio

from embedding import seq_embedding
from AIHDX import prediction

from shiny import App, Inputs, Outputs, Session, render, ui, reactive
from shiny.types import FileInfo

app_ui = ui.page_fluid(
    {"id": "main-content"},
    ui.h2("Predicting protein dynamics with AI-HDX, Try it!"),
    ui.input_text_area("username", "Input your username (i.e. john.smith):"),
    ui.output_ui("show_username"),
    # ui.output_plot("myplot"),
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
    ui.input_action_button("btn", "Submit your data"), # add a click button
    ui.output_table("result"), # output a table
    ui.download_button("downloadData", "Download"), # add a click button for downloading the prediction

)

def server(input: Inputs, output: Outputs, session: Session):
    

    @output
    @render.text
    def show_username():
        username = input.username()
        return "Your username is: " + username
    @output
    @render.ui
    def contents1():
        if input.file1() is None:
            return "Please upload a csv file"
        f: list[FileInfo] = input.file1()
        df = pd.read_csv(f[0]["datapath"])

        # save the data to disk
        username = input.username()
        filename_csv = username + str(time()) + ".csv"
        os.makedirs('AI-HDX_app/save_file', exist_ok=True)
        path = os.path.abspath(os.getcwd())
        df.to_csv(os.path.join(path + "/AI-HDX_app/save_file", filename_csv), header=False, index=False)

        #save file to variable
        df_csv = pd.read_csv(os.path.join(path + "/AI-HDX_app/save_file", filename_csv))

        # return visible table with data
        return None #ui.HTML(df.to_html(classes="table table-striped"))

    @output
    @render.ui

    def contents2():
        if input.file2() is None:
            return "Please upload a txt file"
        f: list[FileInfo] = input.file2()
        df = pd.read_csv(f[0]["datapath"], delimiter="\t")
        
        # save the data to disk
        username = input.username()
        filename_txt = username + str(time()) + ".txt"
        os.makedirs('AI-HDX_app/save_file', exist_ok=True)
        path = os.path.abspath(os.getcwd())
        df.to_csv(os.path.join(path + "/AI-HDX_app/save_file", filename_txt), sep='\t', header=False, index=False)
        
        #save file to variable
        df_txt = pd.read_csv(os.path.join(path + "/AI-HDX_app/save_file", filename_txt), delimiter="\t")

        # return visible table with data
        return None #ui.HTML(df.to_html(classes="table table-striped"))

    # The @reactive.event() causes the function to run only when input.btn is
    # invalidated.
    @reactive.Effect
    @reactive.event(input.btn)
    def _():

        print("Submit your data!")
        

    # This output updates only when input.btn is invalidated.
    @output
    @render.table
    @reactive.event(input.btn)
    def result():
        #run embedding using the input file
        global username     
        username = input.username()  
        path = os.path.abspath(os.getcwd())
        filename_csv = glob.glob(path + "/AI-HDX_app/save_file/"+ username +"*.csv")
        filename_txt = glob.glob(path + "/AI-HDX_app/save_file/"+ username +"*.txt")
        prot= filename_csv[0]
        df = filename_txt[0]
        prot1,df1=seq_embedding(prot, df)
        # make prediction table a global variable
        global output_df
        # run prediction
        output_df = prediction(prot1, df1)
        return output_df
    
    # downloads a file when download button is clicked
    @session.download(        
        filename=lambda: f"AI-HDX_results-{username}-{time()}.csv",

    )

    # writes prediction into the downloaded file
    async def downloadData():       
        yield output_df.to_string(index=False)
        
app = App(app_ui, server)
