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
import sys
if 'D:\\Profil\\Documents\\GitHub\\Fiberphotometry_analysis' not in sys.path:
    sys.path.append('D:\\Profil\\Documents\\GitHub\\Fiberphotometry_analysis')
from Fiberpho_loader import analysis_path, data_path

#%%

#session = 'HC'
mouse = 'C3m'

exp_path = analysis_path / 'Plethysmo'
#session_path = exp_path / f'{session}'
data_path_exp = data_path / '20230128_AliceF_CA2b5Plethysmo'
#mouse_path = Path(f'K:/Alice/Fiber/202110_CA2db2/Analysis/Plethysmo/{session}/{mouse}')

rawdata_path = data_path_exp / f'{mouse}_1.csv'
plethys_df = pd.read_csv(rawdata_path, skiprows=1, usecols=['Time(s)','AIn-4'])
plethys_df = plethys_df.loc[[i for i in range(0,len(plethys_df),600)]] #downsample plethys_df
fiberpho_path = data_path_exp / f'{mouse}_1_dFFfilt.csv'
fiberpho_df = pd.read_csv(fiberpho_path)


if len(fiberpho_df.columns) == 5:
    fiberpho_df.drop(columns='Unnamed: 4', inplace = True)
    fiberpho_df.interpolate(methode = 'nearest', inplace = True)

app = Dash(__name__)

fig = px.line(plethys_df, x='Time(s)', y='AIn-4')
#fig = px.line(fiberpho_df, x='Time(s)', y='Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0-Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass')

app.layout = html.Div([
    html.H4(f'Plethysmo {mouse}'),
    dcc.Graph(
        id='plethysmo',
        figure=fig
    )
])

if __name__ == '__main__':
    app.run_server(debug=False, use_reloader=False)


# #%%to fuse data for A3f
# mouse = 'A3f'
# import numpy as np
# import matplotlib.pyplot as plt

# rawdata1 = pd.read_csv(data_path_exp / f'{mouse}_1old.csv', skiprows=1)
# fiberpho1 = pd.read_csv(data_path_exp / f'{mouse}_1old_dFFfilt.csv')

# rawdata2 = pd.read_csv(data_path_exp / f'{mouse}_2old.csv', skiprows=1)
# fiberpho2 = pd.read_csv(data_path_exp / f'{mouse}_2old_dFFfilt.csv')

# fiberpho2.drop(columns='Unnamed: 4', inplace = True)
# fiberpho2.interpolate(methode = 'nearest', inplace = True)

# rawdata = pd.concat([rawdata1, rawdata2], ignore_index=True)
# fiberpho = pd.concat([fiberpho1, fiberpho2], ignore_index=True)

# #%%
# timevector_raw = np.linspace(0.000083,468.789893+407.408737+0.000083,len(rawdata1)+len(rawdata2))
# timevector_dFF = np.linspace(0.074534,468.574508+407.374624+0.074534,len(fiberpho1)+len(fiberpho2))

# fiberpho['Time(s)'] = timevector_dFF
# rawdata['Time(s)'] = timevector_raw

# plt.plot('Time(s)', 
#          'Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0-Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass', 
#          data = fiberpho1)

# plt.plot('Time(s)', 
#          'Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0-Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass', 
#          data = fiberpho2)

# plt.plot('Time(s)', 
#          'Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0-Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass', 
#          data = fiberpho)

# plt.plot('Time(s)', 
#          'AIn-4', 
#          data = rawdata)

# fiberpho.to_csv(data_path_exp / f'{mouse}_1_dFFfilt.csv')
# #rawdata.to_csv(data_path_exp / f'{mouse}_1.csv')




