import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
from plotly import figure_factory as ff
import matplotlib.pyplot as plt
import io
import plotly.graph_objs as go
import numpy as np
import base64
import FFD_solver
import PFM_solver

image_filename = 'schneider_LIO_White_RGB.png'
encoded_image = base64.b64encode(open(image_filename, 'rb').read())
life_green = '#3DCD58'


def title_header():
    return html.Div(children=[html.Img(src='data:image/png;base64,{}'.format(encoded_image.decode()),
                                       style={'width': 200, 'height': 55.1428571429}),
                              html.Div(style={"flexGrow": 10},
                                       children=[html.H1(
                                           children='FFD/PFM Test-bed',
                                           style={
                                               'textAlign': 'left',
                                               'color': 'white',
                                               'margin': 0,
                                               'padding-left': '30%'
                                           })]
                                       )
                              ],
                    style={'display': 'flex',
                           'backgroundColor': life_green,
                           'padding': '1rem 0.5rem 1rem 0.5rem',
                           'marginBottom': '2rem'})


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
# Loading screen CSS
app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/brPBPO.css"})
# Title of the browser tab
app.title = 'CFD/PFM Test-bed'

app.layout = html.Div(children=[

    title_header(),

    html.Div(children=[
        html.Button(children='Solve',
                    id='button',
                    n_clicks=0,
                    style={'display': 'flex', 'justify-content': 'center',
                           'width': '150px', 'font-size': '18px',
                           'border-radius': '50%'}), ],
        style={'display': 'flex', 'justify-content': 'center'}),

    html.Div(className='twelve columns',
             children=[

                 html.Div(style={'style': 'flex', 'flex-direction': 'column', 'padding-left': '0%'},
                          children=[

                              html.H5(children='Design Parameters',
                                      style={'textAlign': 'left', 'color': 'grey',
                                             'borderBottom': 'solid 1px grey'}),

                              html.Label('Input - 1', style={'font-weight': 'bold', 'textAlign': 'center'}),

                              html.Div(children=[
                                  html.Label(children='Velocity (m/s)',
                                             style={'display': 'inline-block', 'padding-left': '10px',
                                                    'padding-right': '10px'}),
                                  dcc.Input(id='v1-id', value=1, type='number', min=0, max=99, step=0.1,
                                            style={'display': 'inline-block', 'margin-right': '10px',
                                                   'width': '75px', 'text-align': 'center'})],
                                  style={'display': 'flex', 'justify-content': 'space-between',
                                         'padding-right': '20px'}),

                              html.Div(children=[
                                  html.Label(children='Temperature (\u00B0C)',
                                             style={'display': 'inline-block', 'padding-left': '10px',
                                                    'padding-right': '10px'}),
                                  dcc.Input(id='T1-id', value=25, type='number', min=0, max=99, step=0.5,
                                            style={'display': 'inline-block', 'margin-right': '10px',
                                                   'width': '75px', 'text-align': 'center'})],
                                  style={'display': 'flex', 'justify-content': 'space-between',
                                         'padding-right': '20px', 'padding-top': '2%'}),

                              html.Div(children=[
                                  html.Label(children='Angle (\u00B0)',
                                             style={'display': 'inline-block', 'padding-left': '10px',
                                                    'padding-right': '10px'}),
                                  dcc.Input(id='theta1-id', value=90, type='number', min=0, max=180, step=0.5,
                                            style={'display': 'inline-block', 'margin-right': '10px',
                                                   'width': '75px', 'text-align': 'center'})],
                                  style={'display': 'flex', 'justify-content': 'space-between',
                                         'padding-right': '20px', 'padding-top': '2%'}),

                              html.Label('Input - 2', style={'font-weight': 'bold', 'textAlign': 'center'}),

                              html.Div(children=[
                                  html.Label(children='Velocity (m/s)',
                                             style={'display': 'inline-block', 'padding-left': '10px',
                                                    'padding-right': '10px'}),
                                  dcc.Input(id='v2-id', value=1, type='number', min=0, max=99, step=0.1,
                                            style={'display': 'inline-block', 'margin-right': '10px',
                                                   'width': '75px', 'text-align': 'center'})],
                                  style={'display': 'flex', 'justify-content': 'space-between',
                                         'padding-right': '20px'}),

                              html.Div(children=[
                                  html.Label(children='Temperature (\u00B0C)',
                                             style={'display': 'inline-block', 'padding-left': '10px',
                                                    'padding-right': '10px'}),
                                  dcc.Input(id='T2-id', value=25, type='number', min=0, max=99, step=0.5,
                                            style={'display': 'inline-block', 'margin-right': '10px',
                                                   'width': '75px', 'text-align': 'center'})],
                                  style={'display': 'flex', 'justify-content': 'space-between',
                                         'padding-right': '20px', 'padding-top': '2%'}),

                              html.Div(children=[
                                  html.Label(children='Angle (\u00B0)',
                                             style={'display': 'inline-block', 'padding-left': '10px',
                                                    'padding-right': '10px'}),
                                  dcc.Input(id='theta2-id', value=90, type='number', min=0, max=180, step=0.5,
                                            style={'display': 'inline-block', 'margin-right': '10px',
                                                   'width': '75px', 'text-align': 'center'})],
                                  style={'display': 'flex', 'justify-content': 'space-between',
                                         'padding-right': '20px', 'padding-top': '2%'}),

                              html.Br(),

                              html.Div(children=[
                                  html.Label(children='Heat Source',
                                             style={'display': 'inline-block', 'padding-left': '10px',
                                                    'padding-right': '10px'}),
                                  dcc.Input(id='q-id', value=0, type='number', min=0, max=999, step=0.5,
                                            style={'display': 'inline-block', 'margin-right': '10px',
                                                   'width': '75px', 'text-align': 'center'})],
                                  style={'display': 'flex', 'justify-content': 'space-between',
                                         'padding-right': '20px'}),

                              html.Label('Perforated Tile', style={'font-weight': 'bold', 'textAlign': 'center'}),

                              html.Div(children=[
                                  html.Label(children='% Open',
                                             style={'display': 'inline-block', 'padding-left': '10px',
                                                    'padding-right': '10px'}),
                                  dcc.Dropdown(id='b-id', value=0.5,
                                               options=[{'label': '25', 'value': 0.25},
                                                        {'label': '50', 'value': 0.5},
                                                        {'label': '75', 'value': 0.75},
                                                        {'label': '100', 'value': 1.0}],
                                               style={'display': 'inline-block', 'margin-right': '10px',
                                                      'width': '75px', 'text-align': 'center'})],
                                  style={'display': 'flex', 'justify-content': 'space-between',
                                         'padding-right': '20px'}),

                              html.Div(children=[
                                  html.Label(children='Distance from bottom (m)',
                                             style={'display': 'inline-block', 'padding-left': '10px',
                                                    'padding-right': '10px'}),
                                  dcc.Dropdown(id='l-id', value=0.5,
                                               options=[{'label': '0.125', 'value': 0.125},
                                                        {'label': '0.25', 'value': 0.25},
                                                        {'label': '0.375', 'value': 0.375},
                                                        {'label': '0.5', 'value': 0.5},
                                                        {'label': '0.625', 'value': 0.625},
                                                        {'label': '0.75', 'value': 0.75},
                                                        {'label': '0.875', 'value': 0.875}],
                                               style={'display': 'inline-block', 'margin-right': '10px',
                                                      'width': '75px', 'text-align': 'center'})],
                                  style={'display': 'flex', 'justify-content': 'space-between',
                                         'padding-right': '20px'}),
                          ],
                          className='three columns'),

                 html.Div(style={'style': 'flex', 'flex-direction': 'column', 'padding-left': '0%'},
                          children=[

                              html.H5(children='Boundary Conditions',
                                      style={'textAlign': 'left', 'color': 'grey',
                                             'borderBottom': 'solid 1px grey'}),

                              html.Label('Left',
                                         style={'font-weight': 'bold', 'textAlign': 'center'}),

                              html.Div(children=[
                                  html.Label(children='y-Velocity (m/s)',
                                             style={'display': 'inline-block', 'padding-left': '10px',
                                                    'padding-right': '10px'}),
                                  dcc.Input(id='v_l-id', value=0, type='number', min=0, max=99, step=0.1,
                                            style={'display': 'inline-block', 'margin-right': '10px',
                                                   'width': '75px', 'text-align': 'center'})],
                                  style={'display': 'flex', 'justify-content': 'space-between',
                                         'padding-right': '20px'}),

                              html.Div(children=[
                                  html.Label(children='Temperature (\u00B0C)',
                                             style={'display': 'inline-block', 'padding-left': '10px',
                                                    'padding-right': '10px'}),
                                  dcc.Input(id='T_l-id', value=0, type='number', min=0, max=99, step=0.5,
                                            style={'display': 'inline-block', 'margin-right': '10px',
                                                   'width': '75px', 'text-align': 'center'})],
                                  style={'display': 'flex', 'justify-content': 'space-between',
                                         'padding-right': '20px', 'padding-top': '2%'}),

                              html.Div(children=[
                                  html.Label(children='Heat Transfer Coefficient (W/m\u00B2K)',
                                             style={'display': 'inline-block', 'padding-left': '10px',
                                                    'padding-right': '10px'}),
                                  dcc.Input(id='h_l-id', value=1e-15, type='number', min=1e-15, max=99, step=0.5,
                                            style={'display': 'inline-block', 'margin-right': '10px',
                                                   'width': '75px', 'text-align': 'center'})],
                                  style={'display': 'flex', 'justify-content': 'space-between',
                                         'padding-right': '20px', 'padding-top': '2%'}),

                              html.Label('Right',
                                         style={'font-weight': 'bold', 'textAlign': 'center'}),

                              html.Div(children=[
                                  html.Label(children='y-Velocity (m/s)',
                                             style={'display': 'inline-block', 'padding-left': '10px',
                                                    'padding-right': '10px'}),
                                  dcc.Input(id='v_r-id', value=0, type='number', min=0, max=99, step=0.1,
                                            style={'display': 'inline-block', 'margin-right': '10px',
                                                   'width': '75px', 'text-align': 'center'})],
                                  style={'display': 'flex', 'justify-content': 'space-between',
                                         'padding-right': '20px'}),

                              html.Div(children=[
                                  html.Label(children='Temperature (\u00B0C)',
                                             style={'display': 'inline-block', 'padding-left': '10px',
                                                    'padding-right': '10px'}),
                                  dcc.Input(id='T_r-id', value=0, type='number', min=0, max=99, step=0.5,
                                            style={'display': 'inline-block', 'margin-right': '10px',
                                                   'width': '75px', 'text-align': 'center'})],
                                  style={'display': 'flex', 'justify-content': 'space-between',
                                         'padding-right': '20px', 'padding-top': '2%'}),

                              html.Div(children=[
                                  html.Label(children='Heat Transfer Coefficient (W/m\u00B2K)',
                                             style={'display': 'inline-block', 'padding-left': '10px',
                                                    'padding-right': '10px'}),
                                  dcc.Input(id='h_r-id', value=1e-15, type='number', min=1e-15, max=99, step=0.5,
                                            style={'display': 'inline-block', 'margin-right': '10px',
                                                   'width': '75px', 'text-align': 'center'})],
                                  style={'display': 'flex', 'justify-content': 'space-between',
                                         'padding-right': '20px', 'padding-top': '2%'}),

                              html.Label('Top',
                                         style={'font-weight': 'bold', 'textAlign': 'center'}),

                              html.Div(children=[
                                  html.Label(children='x-Velocity (m/s)',
                                             style={'display': 'inline-block', 'padding-left': '10px',
                                                    'padding-right': '10px'}),
                                  dcc.Input(id='u_t-id', value=0, type='number', min=0, max=99, step=0.1,
                                            style={'display': 'inline-block', 'margin-right': '10px',
                                                   'width': '75px', 'text-align': 'center'})],
                                  style={'display': 'flex', 'justify-content': 'space-between',
                                         'padding-right': '20px'}),

                              html.Div(children=[
                                  html.Label(children='Temperature (\u00B0C)',
                                             style={'display': 'inline-block', 'padding-left': '10px',
                                                    'padding-right': '10px'}),
                                  dcc.Input(id='T_t-id', value=0, type='number', min=0, max=99, step=0.5,
                                            style={'display': 'inline-block', 'margin-right': '10px',
                                                   'width': '75px', 'text-align': 'center'})],
                                  style={'display': 'flex', 'justify-content': 'space-between',
                                         'padding-right': '20px', 'padding-top': '2%'}),

                              html.Div(children=[
                                  html.Label(children='Heat Transfer Coefficient (W/m\u00B2K)',
                                             style={'display': 'inline-block', 'padding-left': '10px',
                                                    'padding-right': '10px'}),
                                  dcc.Input(id='h_t-id', value=1e-15, type='number', min=1e-15, max=99, step=0.5,
                                            style={'display': 'inline-block', 'margin-right': '10px',
                                                   'width': '75px', 'text-align': 'center'})],
                                  style={'display': 'flex', 'justify-content': 'space-between',
                                         'padding-right': '20px', 'padding-top': '2%'}),

                              html.Label('Bottom',
                                         style={'font-weight': 'bold', 'textAlign': 'center'}),

                              html.Div(children=[
                                  html.Label(children='x-Velocity (m/s)',
                                             style={'display': 'inline-block', 'padding-left': '10px',
                                                    'padding-right': '10px'}),
                                  dcc.Input(id='u_b-id', value=0, type='number', min=0, max=99, step=0.1,
                                            style={'display': 'inline-block', 'margin-right': '10px',
                                                   'width': '75px', 'text-align': 'center'})],
                                  style={'display': 'flex', 'justify-content': 'space-between',
                                         'padding-right': '20px'}),

                              html.Div(children=[
                                  html.Label(children='Temperature (\u00B0C)',
                                             style={'display': 'inline-block', 'padding-left': '10px',
                                                    'padding-right': '10px'}),
                                  dcc.Input(id='T_b-id', value=0, type='number', min=0, max=99, step=0.5,
                                            style={'display': 'inline-block', 'margin-right': '10px',
                                                   'width': '75px', 'text-align': 'center'})],
                                  style={'display': 'flex', 'justify-content': 'space-between',
                                         'padding-right': '20px', 'padding-top': '2%'}),

                              html.Div(children=[
                                  html.Label(children='Heat Transfer Coefficient (W/m\u00B2K)',
                                             style={'display': 'inline-block', 'padding-left': '10px',
                                                    'padding-right': '10px'}),
                                  dcc.Input(id='h_b-id', value=1e-15, type='number', min=1e-15, max=99, step=0.5,
                                            style={'display': 'inline-block', 'margin-right': '10px',
                                                   'width': '75px', 'text-align': 'center'})],
                                  style={'display': 'flex', 'justify-content': 'space-between',
                                         'padding-right': '20px', 'padding-top': '2%'}),
                          ],
                          className='three columns'),

                 html.Div(style={'style': 'flex', 'flex-direction': 'column', 'padding-left': '0%'},
                          children=[

                              html.H5(children='Modeling - General',
                                      style={'textAlign': 'left', 'color': 'grey',
                                             'borderBottom': 'solid 1px grey'}),

                              html.Div(children=[
                                  html.Label(children='Solver',
                                             style={'display': 'inline-block', 'padding-left': '10px',
                                                    'padding-right': '10px'}),
                                  dcc.Dropdown(id='solver-id', value='CFD',
                                               options=[{'label': 'CFD', 'value': 'CFD'},
                                                        {'label': 'PFM', 'value': 'PFM'}],
                                               style={'display': 'inline-block',
                                                      'margin-right': '10px',
                                                      'width': '75px', 'text-align': 'center'})],
                                  style={'display': 'flex', 'justify-content': 'space-between',
                                         'padding-right': '20px'}),

                              html.Div(children=[
                                  html.Label(children='Grid Cells',
                                             style={'display': 'inline-block', 'padding-left': '10px',
                                                    'padding-right': '10px'}),
                                  dcc.Dropdown(id='N-id', value=16,
                                               options=[{'label': '8', 'value': 8},
                                                        {'label': '16', 'value': 16},
                                                        {'label': '32', 'value': 32},
                                                        {'label': '64', 'value': 64}],
                                               style={'display': 'inline-block',
                                                      'margin-right': '10px',
                                                      'width': '75px', 'text-align': 'center'})],
                                  style={'display': 'flex', 'justify-content': 'space-between',
                                         'padding-right': '20px'}),

                              html.Div(children=[
                                  html.Label(children='SOR Factor',
                                             style={'display': 'inline-block', 'padding-left': '10px',
                                                    'padding-right': '10px'}),
                                  dcc.Input(id='omega-id', value=1.7, type='number', min=0, max=2.0,
                                            step=0.1,
                                            style={'display': 'inline-block', 'margin-right': '10px',
                                                   'width': '75px', 'text-align': 'center'})],
                                  style={'display': 'flex', 'justify-content': 'space-between',
                                         'padding-right': '20px'}),

                              html.Div(children=[
                                  html.Label(children='Temperature ?',
                                             style={'display': 'inline-block', 'padding-left': '10px',
                                                    'padding-right': '10px'}),
                                  dcc.Dropdown(id='temp-id', value='Y',
                                               options=[{'label': 'Y', 'value': 'Y'},
                                                        {'label': 'N', 'value': 'N'}],
                                               style={'display': 'inline-block',
                                                      'margin-right': '10px',
                                                      'width': '75px', 'text-align': 'center'})],
                                  style={'display': 'flex', 'justify-content': 'space-between',
                                         'padding-right': '20px', 'padding-top': '2%'}),

                              html.Div(id='ref_temp',
                                       children=[
                                  html.Label(children='Reference Temp (\u00B0C)',
                                             style={'display': 'inline-block', 'padding-left': '10px',
                                                    'padding-right': '10px'}),
                                  dcc.Input(id='T_ref-id', value=25, type='number', min=1, max=999, step=1,
                                            style={'display': 'inline-block',
                                                   'margin-right': '10px',
                                                   'width': '75px', 'text-align': 'center'})],
                                  style={'display': 'flex', 'justify-content': 'space-between',
                                         'padding-right': '20px'}),

                              html.Div(children=[
                                  html.H5(children='CFD',
                                          style={'textAlign': 'left', 'color': 'grey',
                                                 'borderBottom': 'solid 1px grey'}),

                                  html.Div(children=[
                                      html.Label(children='Time Period (s)',
                                                 style={'display': 'inline-block', 'padding-left': '10px',
                                                        'padding-right': '10px'}),
                                      dcc.Input(id='tmax-id', value=1.0, type='number', min=1.0, max=99, step=1,
                                                style={'display': 'inline-block',
                                                       'margin-right': '10px',
                                                       'width': '75px', 'text-align': 'center'})],
                                      style={'display': 'flex', 'justify-content': 'space-between',
                                             'padding-right': '20px'}),

                                  html.Div(children=[
                                      html.Label(children='Time Step Multiplier',
                                                 style={'display': 'inline-block', 'padding-left': '10px',
                                                        'padding-right': '10px'}),
                                      dcc.Input(id='t_mul-id', value=1.0, type='number', min=0.05, max=50, step=0.05,
                                                style={'display': 'inline-block',
                                                       'margin-right': '10px',
                                                       'width': '75px', 'text-align': 'center'})],
                                      style={'display': 'flex', 'justify-content': 'space-between',
                                             'padding-right': '20px', 'padding-top': '2%'}),

                                  html.Div(children=[
                                      html.Label(children='Velocity Iterations',
                                                 style={'display': 'inline-block', 'padding-left': '10px',
                                                        'padding-right': '10px'}),
                                      dcc.Input(id='vel_iter-id', value=5, type='number', min=1, max=99, step=1,
                                                style={'display': 'inline-block',
                                                       'margin-right': '10px',
                                                       'width': '75px', 'text-align': 'center'})],
                                      style={'display': 'flex', 'justify-content': 'space-between',
                                             'padding-right': '20px', 'padding-top': '2%'}),

                                  html.Div(children=[
                                      html.Label(children='Temp Iterations',
                                                 style={'display': 'inline-block', 'padding-left': '10px',
                                                        'padding-right': '10px'}),
                                      dcc.Input(id='T_iter-id', value=5, type='number', min=1, max=99, step=1,
                                                style={'display': 'inline-block',
                                                       'margin-right': '10px',
                                                       'width': '75px', 'text-align': 'center'})],
                                      style={'display': 'flex', 'justify-content': 'space-between',
                                             'padding-right': '20px', 'padding-top': '2%'}),

                                  html.Div(children=[
                                      html.Label(children='Project Iterations',
                                                 style={'display': 'inline-block', 'padding-left': '10px',
                                                        'padding-right': '10px'}),
                                      dcc.Input(id='P_iter-id', value=50, type='number', min=1, max=999, step=1,
                                                style={'display': 'inline-block',
                                                       'margin-right': '10px',
                                                       'width': '75px', 'text-align': 'center'})],
                                      style={'display': 'flex', 'justify-content': 'space-between',
                                             'padding-right': '20px', 'padding-top': '2%'}),
                              ],
                                  id='CFD_parameters'),

                              # PFM Parameters

                              html.Div(children=[
                                  html.H5(children='PFM',
                                          style={'textAlign': 'left', 'color': 'grey',
                                                 'borderBottom': 'solid 1px grey'}),

                                  html.Div(children=[
                                      html.Label(children='Outer Iterations',
                                                 style={'display': 'inline-block', 'padding-left': '10px',
                                                        'padding-right': '10px'}),
                                      dcc.Input(id='PFM_iter-id', value=5, type='number', min=1, max=99, step=1,
                                                style={'display': 'inline-block',
                                                       'margin-right': '10px',
                                                       'width': '75px', 'text-align': 'center'})],
                                      style={'display': 'flex', 'justify-content': 'space-between',
                                             'padding-right': '20px'}),

                                  html.Div(children=[
                                      html.Label(children='Velocity Potential Iterations',
                                                 style={'display': 'inline-block', 'padding-left': '10px',
                                                        'padding-right': '10px'}),
                                      dcc.Input(id='VP_iter-id', value=75, type='number', min=1, max=999, step=1,
                                                style={'display': 'inline-block',
                                                       'margin-right': '10px',
                                                       'width': '75px', 'text-align': 'center'})],
                                      style={'display': 'flex', 'justify-content': 'space-between',
                                             'padding-right': '20px', 'padding-top': '2%'}),

                                  html.Div(children=[
                                      html.Label(children='Temperature Iterations',
                                                 style={'display': 'inline-block', 'padding-left': '10px',
                                                        'padding-right': '10px'}),
                                      dcc.Input(id='T_iter-id2', value=100, type='number', min=1, max=999, step=1,
                                                style={'display': 'inline-block',
                                                       'margin-right': '10px',
                                                       'width': '75px', 'text-align': 'center'})],
                                      style={'display': 'flex', 'justify-content': 'space-between',
                                             'padding-right': '20px', 'padding-top': '2%'}),
                              ], id='PFM_parameters'),

                              html.Div(children=[
                                  html.H5(children='Monitor Points',
                                          style={'textAlign': 'left', 'color': 'grey',
                                                 'borderBottom': 'solid 1px grey'}),

                                  html.Div(children=[
                                      html.Label(children='Monitor Point x (m)',
                                                 style={'display': 'inline-block', 'padding-left': '10px',
                                                        'padding-right': '10px'}),
                                      dcc.Input(id='mon_x-id', value=0.5, type='number', min=0.0, max=1.0, step=0.1,
                                                style={'display': 'inline-block',
                                                       'margin-right': '10px',
                                                       'width': '75px', 'text-align': 'center'})],
                                      style={'display': 'flex', 'justify-content': 'space-between',
                                             'padding-right': '20px'}),

                                  html.Div(children=[
                                      html.Label(children='Monitor Point y (m)',
                                                 style={'display': 'inline-block', 'padding-left': '10px',
                                                        'padding-right': '10px'}),
                                      dcc.Input(id='mon_y-id', value=0.5, type='number', min=0.0, max=1.0, step=0.1,
                                                style={'display': 'inline-block',
                                                       'margin-right': '10px',
                                                       'width': '75px', 'text-align': 'center'})],
                                      style={'display': 'flex', 'justify-content': 'space-between',
                                             'padding-right': '20px', 'padding-top': '2%'}),
                              ]),

                          ],
                          className='three columns'),

                 html.Div(className='three columns',
                          children=[
                              html.Div(className='twelve columns',
                                       style={'display': 'flex', 'flex-direction': 'column',
                                              'padding-left': '0%'},
                                       children=[
                                           html.H5(children='Results',
                                                   style={'textAlign': 'left', 'color': 'grey',
                                                          'borderBottom': 'solid 1px grey'}),
                                           dcc.Dropdown(id='Q_graph_layout-id',
                                                        # options=[{'label': 'Velocity Vector Plot', 'value': 1},
                                                        #          {'label': 'x-Velocity Contour Plot', 'value': 2},
                                                        #          {'label': 'y-velocity Contour Plot', 'value': 3},
                                                        #          {'label': 'Temperature Contour Plot', 'value': 4},
                                                        #          {'label': 'Pressure Contour Plot', 'value': 5},
                                                        #          {'label': 'Residual Plot', 'value': 6},
                                                        #          {'label': 'Mass Balance Contour Plot', 'value': 7}],
                                                        value=1),
                                           dcc.Graph(id='Q_graph',
                                                     # figure=quiver_fig,
                                                     style={'style': 'flex', 'height': '300px'},
                                                     config={'displayModeBar': False}),
                                           html.Img(id='mpl-id',
                                                    style={'display': 'block', 'height': '100%', 'width': '100%'}),
                                           html.Label('Computational Time (s):'),
                                           html.Label(id='time_complex-id')
                                       ]),
                          ]),
             ]),
])


# Dynamic display - Solver options
@app.callback(
    [Output('CFD_parameters', 'style'), Output('PFM_parameters', 'style')],
    [Input('solver-id', 'value')])
def hide_CFD_PFM(solver):
    if solver == 'CFD':
        return [{'display': 'block'}, {'display': 'none'}]
    else:
        return [{'display': 'none'}, {'display': 'block'}]

# Dynamic display - Reference temperature and temperature contour
@app.callback(
    [Output('ref_temp', 'style'), Output('Q_graph_layout-id', 'options')],
    [Input('temp-id', 'value')])
def hide_ref_temp(solve_temp):
    if solve_temp == 'Y':
        return [{'display': 'flex', 'justify-content': 'space-between', 'padding-right': '20px'},
                [{'label': 'Velocity Vector Plot', 'value': 1},
                 {'label': 'x-Velocity Contour Plot', 'value': 2},
                 {'label': 'y-velocity Contour Plot', 'value': 3},
                 {'label': 'Temperature Contour Plot', 'value': 4},
                 {'label': 'Pressure Contour Plot', 'value': 5},
                 {'label': 'Residual Plot', 'value': 6},
                 {'label': 'Mass Balance Contour Plot', 'value': 7}]]
    else:
        return [{'display': 'none'},
                [{'label': 'Velocity Vector Plot', 'value': 1},
                 {'label': 'x-Velocity Contour Plot', 'value': 2},
                 {'label': 'y-velocity Contour Plot', 'value': 3},
                 {'label': 'Pressure Contour Plot', 'value': 5},
                 {'label': 'Residual Plot', 'value': 6},
                 {'label': 'Mass Balance Contour Plot', 'value': 7}]]

# Calculation
@app.callback(
    [Output('Q_graph', 'figure'), Output('mpl-id', 'src'), Output('time_complex-id', 'children')],
    [Input('button', 'n_clicks')],
    [State('v1-id', 'value'), State('T1-id', 'value'), State('theta1-id', 'value'),
     State('v2-id', 'value'), State('T2-id', 'value'), State('theta2-id', 'value'),
     State('q-id', 'value'), State('b-id', 'value'), State('l-id', 'value'),
     State('v_l-id', 'value'), State('T_l-id', 'value'), State('h_l-id', 'value'),
     State('v_r-id', 'value'), State('T_r-id', 'value'), State('h_r-id', 'value'),
     State('u_t-id', 'value'), State('T_t-id', 'value'), State('h_t-id', 'value'),
     State('u_b-id', 'value'), State('T_b-id', 'value'), State('h_b-id', 'value'),
     State('solver-id', 'value'), State('N-id', 'value'), State('omega-id', 'value'), State('temp-id', 'value'),
     State('tmax-id', 'value'), State('t_mul-id', 'value'),
     State('vel_iter-id', 'value'), State('T_iter-id', 'value'), State('P_iter-id', 'value'),
     State('mon_x-id', 'value'), State('mon_y-id', 'value'), State('T_ref-id', 'value'),
     State('PFM_iter-id', 'value'), State('VP_iter-id', 'value'), State('T_iter-id2', 'value'),
     State('Q_graph_layout-id', 'value')
     ])
def cfd_results(n_clicks, v1, T1, theta1, v2, T2, theta2, q_total, b, L, v_l, T_l, h_l, v_r, T_r, h_r, u_t, T_t, h_t, u_b, T_b,
                h_b, solver, N, omega, solve_temp, tmax, tmul, vel_iter, T_iter, P_iter, mon_x, mon_y, T_ref, PFM_iter,
                VP_iter, T_iter_2, Q_graph_layout):
    design_inputs = []
    solver_config = []

    input_1 = []
    input_1.append(v1)
    input_1.append(T1)
    input_1.append(theta1)
    design_inputs.append(input_1)

    input_2 = []
    input_2.append(v2)
    input_2.append(T2)
    input_2.append(theta2)
    design_inputs.append(input_2)

    design_inputs.append(q_total)

    opening = []
    opening.append(b)
    opening.append(L)
    design_inputs.append(opening)

    left_wall = []
    left_wall.append(v_l)
    left_wall.append(T_l)
    left_wall.append(h_l)
    design_inputs.append(left_wall)

    right_wall = []
    right_wall.append(v_r)
    right_wall.append(T_r)
    right_wall.append(h_r)
    design_inputs.append(right_wall)

    top_wall = []
    top_wall.append(u_t)
    top_wall.append(T_t)
    top_wall.append(h_t)
    design_inputs.append(top_wall)

    bottom_wall = []
    bottom_wall.append(u_b)
    bottom_wall.append(T_b)
    bottom_wall.append(h_b)
    design_inputs.append(bottom_wall)

    solver_config.append(N)
    solver_config.append(omega)

    temp = []
    temp.append(solve_temp)
    temp.append(T_ref)
    solver_config.append(temp)

    CFD_parameter = []
    CFD_parameter.append(tmax)
    CFD_parameter.append(tmul)
    CFD_parameter.append(vel_iter)
    CFD_parameter.append(T_iter)
    CFD_parameter.append(P_iter)
    solver_config.append(CFD_parameter)

    monitor = []
    monitor.append(mon_x)
    monitor.append(mon_y)
    solver_config.append(monitor)

    PFM_parameter = []
    PFM_parameter.append(PFM_iter)
    PFM_parameter.append(VP_iter)
    PFM_parameter.append(T_iter_2)
    solver_config.append(PFM_parameter)

    # Time complexity calculation
    if solver == 'CFD':
        CFDSolve = FFD_solver.CFD(design_inputs, solver_config)
        U_col, V_col, T_col, P_col, timestamp, monitor_data, mass = CFDSolve.Solve_CFD()
    else:
        PFMSolve = PFM_solver.PFM(design_inputs, solver_config)
        U_col, V_col, T_col, P_col, timestamp, monitor_data, mass = PFMSolve.Solve_PFM()

    U_col = np.round_(U_col, 4)
    V_col = np.round_(V_col, 4)
    T_col = np.round_(T_col, 4)
    P_col = np.round_(P_col, 4)
    # mass = np.round_(mass, 4)

    x, y = np.meshgrid(np.arange(1, N + 1), np.arange(1, N + 1))

    # Vector plot ---------------------------------------
    fig1, ax1 = plt.subplots()
    M = np.hypot(U_col, V_col)
    Q = ax1.quiver(x, y, U_col, V_col, M)

    vect_img = io.BytesIO()
    plt.savefig(vect_img, format='png', bbox_inches='tight')
    vect_img.seek(0)
    vect_fig = base64.b64encode(vect_img.getvalue())

    vect_src='data:image/png;base64,{}'.format(vect_fig.decode('utf-8'))

    # Contour plots --------------------------------------
    # U_velocity plot
    plt.figure()
    u_vel = plt.contourf(x, y, V_col, 50)
    plt.colorbar(u_vel)

    u_img = io.BytesIO()
    plt.savefig(u_img, format='png', bbox_inches='tight')
    u_img.seek(0)
    u_fig = base64.b64encode(u_img.getvalue())

    u_src='data:image/png;base64,{}'.format(u_fig.decode('utf-8'))

    # V_velocity plot
    plt.figure()
    v_vel = plt.contourf(x, y, V_col, 50)
    plt.colorbar(v_vel)

    v_img = io.BytesIO()
    plt.savefig(v_img, format='png', bbox_inches='tight')
    v_img.seek(0)
    v_fig = base64.b64encode(v_img.getvalue())

    v_src='data:image/png;base64,{}'.format(v_fig.decode('utf-8'))

    # Pressure plot
    plt.figure()
    p = plt.contourf(x, y, P_col, 50)
    plt.colorbar(p)

    p_img = io.BytesIO()
    plt.savefig(p_img, format='png', bbox_inches='tight')
    p_img.seek(0)
    p_fig = base64.b64encode(p_img.getvalue())

    p_src='data:image/png;base64,{}'.format(p_fig.decode('utf-8'))

    # Temperature plot
    plt.figure()
    temp = plt.contourf(x, y, T_col, 100)
    plt.colorbar(temp)

    temp_img = io.BytesIO()
    plt.savefig(temp_img, format='png', bbox_inches='tight')
    temp_img.seek(0)
    temp_fig = base64.b64encode(temp_img.getvalue())

    temp_src='data:image/png;base64,{}'.format(temp_fig.decode('utf-8'))

    # Residual plot
    plt.figure()
    plt.plot(monitor_data[0], monitor_data[1], label='x-velocity')
    plt.plot(monitor_data[0], monitor_data[2], label='y-velocity')
    plt.plot(monitor_data[0], monitor_data[3], label='Pressure')
    plt.plot(monitor_data[5], monitor_data[4], label='Temp')
    plt.legend()

    mon_img = io.BytesIO()
    plt.savefig(mon_img, format='png', bbox_inches='tight')
    mon_img.seek(0)
    mon_fig = base64.b64encode(mon_img.getvalue())

    mon_src='data:image/png;base64,{}'.format(mon_fig.decode('utf-8'))

    # Mass balance plot
    plt.figure()
    mass_balance = plt.contourf(x, y, mass, 100)
    plt.colorbar(mass_balance)

    mass_img = io.BytesIO()
    plt.savefig(mass_img, format='png', bbox_inches='tight')
    mass_img.seek(0)
    mass_fig = base64.b64encode(mass_img.getvalue())

    mass_src='data:image/png;base64,{}'.format(mass_fig.decode('utf-8'))

    if Q_graph_layout == 1:
        # Velocity Vector Graph
        quiver_fig = ff.create_quiver(x, y, U_col, V_col,
                                      scale=1,
                                      arrow_scale=0.4,
                                      name='quiver',
                                      line=dict(width=1))
        quiver_fig.layout.margin.update({'l': 25, 't': 10, 'r': 10, 'b': 20})
        return [quiver_fig, vect_src, timestamp]

    elif Q_graph_layout == 2:
        # U velocity contour plot
        U_contour = {'data': [go.Contour(z=U_col, x=np.arange(1, N + 1), y=np.arange(1, N + 1), line_smoothing=1.3, ncontours=50)],
                     'layout': {'margin': {'l': 25, 't': 50, 'r': 0, 'b': 30}}}
        return [U_contour, u_src, timestamp]

    elif Q_graph_layout == 3:
        # V velocity contour plot
        V_contour = {'data': [go.Contour(z=V_col, x=np.arange(1, N + 1), y=np.arange(1, N + 1), line_smoothing=1.3, ncontours=50)],
                     'layout': {'margin': {'l': 25, 't': 50, 'r': 0, 'b': 30}}}
        return [V_contour, v_src, timestamp]

    elif Q_graph_layout == 4:
        # Temperature contour plot
        T_contour = {'data': [go.Contour(z=T_col, x=np.arange(1, N + 1), y=np.arange(1, N + 1), line_smoothing=1.3, ncontours=50)],
                     'layout': {'margin': {'l': 25, 't': 50, 'r': 0, 'b': 30}}}
        return [T_contour, temp_src, timestamp]

    elif Q_graph_layout == 5:
        # Pressure contour plot
        P_contour = {'data': [go.Contour(z=P_col, x=np.arange(1, N + 1), y=np.arange(1, N + 1), line_smoothing=1.3, ncontours=50)],
                     'layout': {'margin': {'l': 20, 't': 50, 'r': 0, 'b': 30}}}
        return [P_contour, p_src, timestamp]

    elif Q_graph_layout == 6:
        # Residual graph
        Res_plot = {'data': [go.Scatter(x=monitor_data[0], y=monitor_data[1], mode='lines'),
                             go.Scatter(x=monitor_data[0], y=monitor_data[2], mode='lines'),
                             go.Scatter(x=monitor_data[0], y=monitor_data[3], mode='lines'),
                             go.Scatter(x=monitor_data[0], y=monitor_data[4], mode='lines')],
                    'layout': {'margin': {'l': 20, 't': 50, 'r': 0, 'b': 30},
                               'legend': {'x': 0, 'y': 1.1, 'orientation': 'h'}}}
        return [Res_plot, mon_src, timestamp]

    elif Q_graph_layout == 7:
        # Mass Balance contour plot
        mass_contour = {'data': [go.Contour(z=mass, x=np.arange(1, N + 1), y=np.arange(1, N + 1), line_smoothing=1.3, ncontours=50)],
                     'layout': {'margin': {'l': 20, 't': 50, 'r': 0, 'b': 30}}}
        return [mass_contour, mass_src, timestamp]


if __name__ == '__main__':
    app.run_server(debug=True, port=8050)
