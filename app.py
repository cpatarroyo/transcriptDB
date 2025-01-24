from shiny import App, Inputs, Outputs, Session, reactive, render, ui
import pandas as pd

app_ui = ui.page_fluid(
    ui.panel_title("Transcriptoma de Scedosporium apiospermum"),
    ui.layout_sidebar(
        ui.sidebar(
            ui.input_text('searchTerm', label='Termino de busqueda'),
            ui.input_select('searchType', label='Buscar por:', choices=['Nombre', 'Anotación', 'Numero EC', 'Categoría Pfam', 'Reacción KEGG']),
            ui.input_switch('filter', label='Significancia'),
            ui.input_numeric('threshold', label='Umbral', min=0, max=1, step=0.01, value=0.05),
            ui.input_action_button('search', label='Buscar')
        ),
        ui.output_data_frame('resultsTab')
    )
)

def server(input: Inputs, output: Outputs, session: Session):
    @reactive.event(input.search)
    def getResults():
        test = pd.DataFrame({'col1': [1, 2], 'col2': [3, 4]})
        return test

    @render.data_frame
    def resultsTab():
        return render.DataGrid(getResults(),selection_mode='rows')
        #return getResults()

app = App(app_ui, server)