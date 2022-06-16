# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 11:24:30 2022

@author: afermigier
"""
# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

from dash import Dash, dcc, html
import plotly.express as px
import pandas as pd
from pathlib import Path
import sys
if 'D:\\Profil\\Documents\\GitHub\\Fiberphotometry_analysis' not in sys.path:
    sys.path.append('D:\\Profil\\Documents\\GitHub\\Fiberphotometry_analysis')
from Fiberpho_loader import analysis_path, data_path

session = 'Novel'
mouse = 'CDm1'

exp_path = analysis_path / 'Plethysmo'
session_path = exp_path / f'{session}'
data_path_exp = data_path / '20211022_AliceF_CA2b2plethysmoNovel'
mouse_path = Path(f'K:/Alice/Fiber/202110_CA2db2/Analysis/Plethysmo/{session}/{mouse}')

rawdata_path = data_path_exp / f'{mouse}_0.csv'
plethys_df = pd.read_csv(rawdata_path, skiprows=1, usecols=['Time(s)','AIn-4'])
plethys_df = plethys_df.loc[[i for i in range(0,len(plethys_df),600)]] #downsample plethys_df
fiberpho_path = data_path_exp / f'{mouse}_0_dFFfilt.csv'
fiberpho_df = pd.read_csv(fiberpho_path)


if len(fiberpho_df.columns) == 5:
    fiberpho_df.drop(columns='Unnamed: 4', inplace = True)
    fiberpho_df.interpolate(methode = 'nearest', inplace = True)

app = Dash(__name__)

fig = px.line(plethys_df, x='Time(s)', y='AIn-4')
#fig = px.line(fiberpho_df, x='Time(s)', y='Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0-Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass')

app.layout = html.Div([
    html.H4(f'Plethysmo {session} {mouse}'),
    dcc.Graph(
        id='plethysmo',
        figure=fig
    )
])

if __name__ == '__main__':
    app.run_server(debug=True)





