import dash
import dash_core_components as dcc
from dash.dependencies import Input, Output, State
import dash_html_components as html
import plotly.graph_objs as go
import scanpy as sc
import pandas as pd
import dash_table
import base64
import uuid
from flask import Flask, send_from_directory
import os
from urllib.parse import quote as urlquote

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

path = 'C:\\code\\bioinformatics\\scRNA-workshop\\scRNA-python-workshop-master\\content'

UPLOAD_DIRECTORY = "C:/code/bioinformatics/rna-seq-static"

if not os.path.exists(UPLOAD_DIRECTORY):
    os.makedirs(UPLOAD_DIRECTORY)

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
server = Flask(__name__)
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
session_id = str(uuid.uuid4())
app.layout = html.Div(children=[

    dcc.Tabs([
        #DATA UPLOAD TAB
        dcc.Tab(label='Data Upload', children=[
            html.Div([html.H2('Data Upload')]),
            dcc.Upload(
                id='upload-data',
                children=html.Div([
                    'Drag and Drop or ',
                    html.A('Select Files')
                ]),
                style={
                    'width': '100%',
                    'height': '60px',
                    'lineHeight': '60px',
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'margin': '10px'
                },
                multiple=True
            ),
            html.Div(id='table-head'),
            html.Div(children=[html.Div(session_id, id='session-id', hidden='hidden')])

        ]),

        dcc.Tab(label='Data Visualization', children=[
            html.Div([html.H2('Data Visualization')]),
            dcc.Dropdown(id='analysis-type',
                         options=[{'label': s, 'value': s} for s in algorithms.keys()],
                         value=[],
                         multi=True),

            dcc.Dropdown(id='meta-type',
                         options=[{'label': m, 'value': m} for m in adata.obs],
                         value=[]),

            html.Div(children=html.Div(id='graphs', className='row'))
        ]),

    ])


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
                                                       showlegend=True,
                                                       legend=dict(
                                                        yanchor="top",
                                                        y=0.99,
                                                        xanchor="left",
                                                        x=0.01
                                                    ))}
        )))

    return graphs



import io
from datetime import datetime
static_path = 'C:/code/bioinformatics/rna-seq-static/'

def parse_file(contents, filename): #returns a table of the inputted data
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        if 'csv' in filename:
            # Assume that the user uploaded a CSV file
            df = pd.read_csv(
                io.StringIO(decoded.decode('utf-8')))
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            df = pd.read_excel(io.BytesIO(decoded))
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])

    df.to_json(static_path + filename.split('.')[0] + '.json')

    return html.Div([
        html.H5(filename),
        html.H6(datetime.now()),

        dash_table.DataTable(
            data=df.head(25).to_dict('records'),
            columns=[{'name': i, 'id': i} for i in df.columns]
        ),

        html.Hr(),  # horizontal line

        # For debugging, display the raw contents provided by the web browser
        html.Div('Raw Content'),
        html.Pre(contents[0:200] + '...', style={
            'whiteSpace': 'pre-wrap',
            'wordBreak': 'break-all'
        })
    ])


@app.callback(Output('table-head', 'children'),
              [Input('upload-data', 'contents')],
              [State('upload-data', 'filename')])
def update_output(list_of_contents, list_of_names):
    if list_of_contents is not None:
        children = [
            parse_file(c, n) for c, n in
            zip(list_of_contents, list_of_names)]
        return children



# GOAL: Update graph's colors without reloading data
# @app.callback(
#     Output('graphs', 'children'), [Input('meta-type', 'value')]
# )
# def apply_meta(meta_type):
#     for graph in


if __name__ == '__main__':
    app.run_server(debug=True)
