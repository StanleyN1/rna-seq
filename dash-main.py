import dash
import dash_core_components as dcc
from dash.dependencies import Input, Output
import dash_html_components as html
import plotly.graph_objs as go
import scanpy as sc
import pandas as pd

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

path = 'C:\\code\\bioinformatics\\scRNA-workshop\\scRNA-python-workshop-master\\content'

adata = sc.read(path + '/data/brain_analysis.h5ad')

pca = adata.obsm['X_pca']
umap = adata.obsm['X_umap']
tsne = adata.obsm['X_tsne']

algorithms = {
    'PCA' : pca,
    'UMAP': umap,
    'TSNE': tsne
}
#COLOR SCHEME, LEGEND,
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div(children=[

    html.Div([html.H2('Data Visualizer')]),
    dcc.Dropdown(id='analysis-type',
                 options=[{'label': s, 'value': s} for s in algorithms.keys()],
                 value=[],
                 multi=True),

    dcc.Dropdown(id='meta-type',
                 options=[{'label': m, 'value': m} for m in adata.obs],
                 value=[]),

    html.Div(children=html.Div(id='graphs', className='row'))

], className='container', style={'width' : '98%', 'margin-left': 10, 'margin-right': 10})

from random import randint
@app.callback(
    Output('graphs', 'children'), [Input('analysis-type', 'value'), Input('meta-type', 'value')]
)
def update_graph(analysis_types, meta_type):
    graphs = []
    meta_colors = {}
    meta_color = []
    if meta_type:
        meta_colors = { #generates random hex color values for the meta data types
            str(m) : "#{:06x}".format(randint(0x3FFFFF, 0xBFFFFF)) for m in adata.obs[str(meta_type)].unique()
        }
        meta_color = [meta_colors[c] for c in list(adata.obs[str(meta_type)])] #converts meta data value to hex value

    for analysis in analysis_types:
        x=[data[0] for data in algorithms[analysis]]
        y=[data[1] for data in algorithms[analysis]]

        data = go.Scatter(
            x=x,
            y=y,
            mode='markers',
            marker=dict(color=meta_color)
        )

        graphs.append(html.Div(dcc.Graph(
            id=analysis,
            figure={'data':[data], 'layout': go.Layout(xaxis=dict(range=[min(x), max(x)]),
                                                       yaxis=dict(range=[min(y), max(y)]),
                                                       margin={'l':50, 'r':1, 't':45, 'b':1},
                                                       title=analysis,
                                                       showlegend=True)}
        )))

    return graphs



# GOAL: Update graph's colors without reloading data
# @app.callback(
#     Output('graphs', 'children'), [Input('meta-type', 'value')]
# )
# def apply_meta(meta_type):
#     for graph in


if __name__ == '__main__':
    app.run_server(debug=True)
