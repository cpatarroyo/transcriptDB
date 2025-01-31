from shiny import App, Inputs, Outputs, Session, reactive, render, ui
import pandas as pd
import matplotlib.pyplot as mpl
import numpy as np
from dbFunc import nameSearch, ecSearch, pfamSearch, keggReactSearch, descriptionSearch

app_ui = ui.page_fluid(
    ui.panel_title("Transcriptoma de Scedosporium apiospermum"),
    ui.layout_sidebar(
        ui.sidebar(
            ui.input_text('searchTerm', label='Termino de busqueda'),
            ui.input_select('searchType', label='Buscar por:', choices=['Nombre', 'Descripción', 'Numero EC', 'Categoría Pfam', 'Reacción KEGG']),
            ui.input_switch('filter', label='Significancia'),
            ui.input_numeric('threshold', label='Umbral', min=0, max=1, step=0.01, value=0.05),
            ui.input_action_button('search', label='Buscar'),
            ui.output_plot('difexp')
        ),
        ui.output_data_frame('resultsTab')
    )
)

def server(input: Inputs, output: Outputs, session: Session):
    @reactive.event(input.search)
    def getResults():
        if input.searchType() == 'Nombre':
            results = nameSearch(input.searchTerm(), signif = input.filter(), signifThr = float(input.threshold()))
        elif input.searchType() == 'Descripción':
            results = descriptionSearch(input.searchTerm(), signif = input.filter(), signifThr = float(input.threshold()))
        elif input.searchType() == 'Numero EC':
            results = ecSearch(input.searchTerm(), signif = input.filter(), signifThr = float(input.threshold()))
        elif input.searchType() == 'Categoría Pfam':
            results = pfamSearch(input.searchTerm(), signif = input.filter(), signifThr = float(input.threshold()))
        elif input.searchType() == 'Reacción KEGG':
            results = keggReactSearch(input.searchTerm(), signif = input.filter(), signifThr = float(input.threshold()))
        return results

    @render.data_frame
    def resultsTab():
        return render.DataGrid(getResults(),selection_mode='rows')

    @render.plot
    def difexp():
        recount = pd.cut(getResults()['LogFC'],[-np.inf,-1,1,np.inf], labels=['RE','NON','OE']).value_counts()
        fig, ax = mpl.subplots()
        ax.pie([recount['RE'], recount['NON'], recount['OE']], labels = [f'RE({recount["RE"]})', f'No({recount["NON"]})', f'SE({recount["OE"]})'], autopct = "%1.0f%%", colors = ['red', 'gray','green'])

app = App(app_ui, server)