# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 15:08:44 2023

To visualize data in Dash using plotly

@author: alice fermigier
"""

#%%
##########
#IMPORTED#
##########

from dash import Dash, dcc, html
import plotly.express as px

#%%
###################
#DEFINED FUNCTIONS#
###################

def visualize(data_df,data_type,exp,session,mouse):
     
    app = Dash(__name__)
    
    if data_type == 'fiberpho':
        fig = px.line(data_df[50:], x='Time(s)', y='470 Deinterleaved')
        #fig = px.line(deinterleaved_df, x='Time(s)', y='405 Deinterleaved')
        
    elif data_type == 'plethysmo':
        data_df = data_df.loc[[i for i in range(0,len(data_df),600)]] #downsample plethysmo data
        fig = px.line(data_df, x='Time(s)', y='AIn-4')
    
    app.layout = html.Div([
        html.H4(f'{exp} {session} {mouse}'),
        dcc.Graph(
            id=f'{exp}',
            figure=fig
        )
    ])
    
    if __name__ == '__main__':
        app.run_server(debug=False, use_reloader=False)
