from dash import Dash, dcc, html, Input, Output, State
import plotly.express as px
import pandas as pd
import json
import numpy as np
import importlib

import modules.common.nomenclature as nom
importlib.reload(nom)
import modules.common.preprocess as pp
importlib.reload(pp)
import modules.behaviour.fear_conditioning as fc
importlib.reload(fc)

from scripts.loader import experiment_path, data_path, proto_df, TIME_BEGIN, batches

#------------------#
exp = 'FCConditioning'
sheet = 'Conditioning'
mouse = '1009'
batch = 4
#------------------# 

datapath_exp_dict = nom.get_experiment_data_path(batches, proto_df, data_path, exp)
data_path_exp = datapath_exp_dict[batch]
pp_path = data_path_exp / 'Preprocessing'
behav_path_exp = data_path_exp / 'Behaviour'
fiberpho_path = pp_path / f'{mouse}_dFF_corrected_final.csv'
fiberpho_df = pd.read_csv(fiberpho_path)
protocol_file = experiment_path / "fear_protocol.xlsx"
proto = fc.parse_protocol_sheet(protocol_file, sheet)

# Path where results will be stored
score_path = behav_path_exp / f'{batch}_{mouse}_manual_shock_scoring.json'

shock_imetronic = np.array([onset for onset, _ in proto["Shock"]])
n_shocks = len(shock_imetronic)
time = fiberpho_df['Time(s)'].values

# Create the Dash app
app = Dash(__name__)

# Create the figure
fig = px.line(fiberpho_df[TIME_BEGIN:], x='Time(s)', y='dFF')

# App layout
app.layout = html.Div([
    html.H4(f'{exp} {mouse}'),
    
    dcc.Graph(
        id='plot',
        figure=fig,
        config={'displayModeBar': True}
    ),
    
    html.Div(id='scorer-message', style={'color': 'black', 'fontWeight': 'bold'}),
    
    html.Button("Save Excluded Regions", id="save-button", n_clicks=0),
    
    dcc.Store(id='scorer-storage', data=[]),
    dcc.Store(id='click-tracker', data=None)
])

# ── Dash app ──────────────────────────────────────────────────────────────
app.layout = html.Div([
    html.H4(f'Shock scorer — Batch {batch}, Mouse {mouse}',
            style={'marginBottom': '4px'}),
    html.P(
        f'Expected Imetronic shock times (s): {np.round(shock_imetronic, 1).tolist()}',
        style={'fontSize': '12px', 'color': '#555'}
    ),
    dcc.Graph(
        id='fiber-plot', figure=fig,
        config={'displayModeBar': True},
        style={'height': '360px'}
    ),
    html.Div(id='score-message',
                style={'color': '#333', 'fontWeight': 'bold',
                    'marginTop': '6px', 'fontSize': '13px'}),
    html.Button('Save & Compute Mapping',
                id='save-btn', n_clicks=0,
                style={'marginTop': '8px'}),
    html.Div(id='save-message',
                style={'color': 'green', 'marginTop': '6px'}),

    dcc.Store(id='clicks-store', data=[]),   # list of clicked Doric times
])

# ── Callback: record click ─────────────────────
@app.callback(
    [Output('score-message','children'),
        Output('clicks-store', 'data')],
    Input('fiber-plot',  'clickData'),
    [State('clicks-store',  'data')]
)
def record_click(click_data, scored_times):
    if not click_data:
        remaining = n_shocks - len(scored_times)
        return f'Click shock onset {len(scored_times)+1} / {n_shocks}', scored_times

    t_clicked = click_data['points'][0]['x']

    # Ignore extra clicks beyond n_shocks
    if len(scored_times) >= n_shocks:
        return f'All {n_shocks} shocks scored. Hit Save (or click to redo).', scored_times

    scored_times = scored_times + [t_clicked]   # immutable append

    if len(scored_times) < n_shocks:
        msg = f'Shock {len(scored_times)} marked at {t_clicked:.3f} s — click shock {len(scored_times)+1}'
    else:
        msg = f'All {n_shocks} shocks marked — check then hit Save.'

    return msg, scored_times

# ── Callback: save + compute ──────────────────────────────────────────
@app.callback(
    Output('save-message', 'children'),
    Input('save-btn',      'n_clicks'),
    State('clicks-store',  'data')
)
def save_scoring(n_clicks, scored_times):
    if n_clicks == 0:
        return ''
    if len(scored_times) != n_shocks:
        return f'[!] Need {n_shocks} clicks, have {len(scored_times)}. Keep clicking.'

    doric_onsets = np.array(scored_times)
    effective_slope, protocol_start = np.polyfit(shock_imetronic, doric_onsets, 1)

    predicted     = protocol_start + effective_slope * shock_imetronic
    residual_rms  = np.sqrt(np.mean((doric_onsets - predicted)**2)) * 1000

    result = {
        'protocol_start':   float(protocol_start),
        'effective_slope':  float(effective_slope),
        'doric_onsets':     doric_onsets.tolist(),
        'imetronic_onsets': shock_imetronic.tolist(),
        'residual_rms_ms':  float(residual_rms),
    }
    with open(score_path, 'w') as f:
        json.dump(result, f, indent=2)

    print(f"[Manual scoring] Saved to {score_path.name}")
    print(f"  protocol_start  = {protocol_start:.4f} s")
    print(f"  effective_slope = {effective_slope:.8f}")
    print(f"  Residual RMS    = {residual_rms:.1f} ms")

    return (f'Saved. protocol_start={protocol_start:.3f} s | '
            f'slope={effective_slope:.7f} | '
            f'residual RMS={residual_rms:.1f} ms')

app.run(debug=False, use_reloader=False)