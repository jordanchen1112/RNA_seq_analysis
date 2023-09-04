import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State
from Interactive_plot import InteractivePlotter
import dash.exceptions

app = dash.Dash(__name__)

# Initial default values for the thresholds
default_p = 0.05
default_fc = 1.5

plotter = InteractivePlotter()  # Assuming you have the InteractivePlotter class defined as above

app.layout = html.Div([
    # Input for gene of interest
    html.Label('Enter query of interest:'),
    dcc.Input(id='interest-input', type='text', value=''),

    # Initial number of clicks is 0
    html.Button('Update Plot', id='submit-button', n_clicks=0),  

    # Sliders for adjusting thresholds
    html.Label('Adjust p-value threshold:'),
    dcc.Slider(
        id='p-slider',
        min=0.01,
        max=0.1,
        step=0.01,
        value=default_p,
        marks={i/100: str(i/100) for i in range(1, 11)},
    ),
    html.Label('Adjust fold change threshold:'),
    dcc.Slider(
        id='fc-slider',
        min=1,
        max=4,
        step=0.1,
        value=default_fc,
        marks={i: str(i) for i in range(1, 5)},
    ),
    # Div to display the volcano plot
    html.Div(dcc.Graph(id='volcano-plot'), style={'width': '1200px', 'height': '1200px'}),
])

@app.callback(
    Output('volcano-plot', 'figure'),
    [Input('submit-button', 'n_clicks')],  # Only the button as an Input
    [State('p-slider', 'value'),  # p-slider and fc-slider as States
     State('fc-slider', 'value'),
     State('interest-input', 'value')]
)
def update_figure(n_clicks, p_value, fold_change, interest):
    if not interest or interest == '':
        interest = None
    fig = plotter.interactive_volcano(p=p_value, fc=fold_change, interest=interest)
    return fig
    
if __name__ == '__main__':
    app.run_server(debug=True, use_reloader=False, port=8050)
