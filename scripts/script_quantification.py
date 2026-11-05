#%%IMPORTED
###########

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import os
import importlib

#import functions
import modules.common.preprocess as pp
importlib.reload(pp)
import modules.common.genplot as gp
importlib.reload(gp)
import modules.common.behavplot as bp
importlib.reload(bp)
import modules.common.statcalc as sc
importlib.reload(sc)
import modules.common.transients as tr
importlib.reload(tr)
import modules.common.nomenclature as nom
importlib.reload(nom)
import modules.behaviour.mouse_position as mp
importlib.reload(mp)
import modules.behaviour.epm as epm
importlib.reload(epm)
import modules.behaviour.camera_processing as cp
importlib.reload(cp)
import modules.common.clean_signal as cs
importlib.reload(cs)
import modules.common.quantification as quantif
importlib.reload(quantif)

from scripts.loader import analysis_path, data_path, proto_df, subjects_df, batches, experiment_path

#filter characteristics
ORDER = 2
CUT_FREQ = 4 #in Hz
#threshold to fuse behaviour if bouts are too close, in secs
THRESH_S = 0
#threshold for PETH : if events are too short do not plot them and do not include them in PETH, in seconds
EVENT_TIME_THRESHOLD = 0

exp = 'Plethysmo'
# if plethysmo, chose concatenated data or not
concat = True

exp_path = analysis_path / exp
data_path_exp = data_path / proto_df.loc[proto_df['Task']==exp, 'Data_path'].values[0]
pp_path = data_path_exp / 'Preprocessing'

# %% Extract transients (count, frequency, mean amplitude, mean width)

repo_path = exp_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
subjects_df['Group'] = subjects_df['Group'].fillna('')

plot_dir = repo_path / 'transient_plots'
plot_dir.mkdir(exist_ok=True)

excluded_subjects_df = pd.read_excel(experiment_path / 'subjects.xlsx',
                                     sheet_name=f'Excluded_{exp}')
excluded_ids = set(excluded_subjects_df['Subject'].astype(str))

# Two detection strategies, see note below
configs = {
    'mad3': dict(threshold_mode='mad', mad_k=3.0),
    'abs':  dict(threshold_mode='absolute', abs_prominence=3.0),   # <-- set in dFF units
}
channels = ['dFF', '560 dFF']

results = {(ch, cfg): {'summary': [], 'events': []} for ch in channels for cfg in configs}

for mouse, batch, group in zip(subjects_df['Subject'], subjects_df['Batch'], subjects_df['Group']):
    print(f"--- {mouse} {batch} {group} ---")
    if str(mouse) in excluded_ids:
        print(f"Mouse {mouse} excluded")
        continue
    if concat:
        fiberbehav_file = repo_path / f'{batch}_{mouse}_fiberbehavconcatnotderived.csv'
    else:
        fiberbehav_file = repo_path / f'{batch}_{mouse}_fiberbehavnotderived.csv'
    if not fiberbehav_file.exists():
        print("  File not found, skipping.")
        continue

    fiberbehav_df = pd.read_csv(fiberbehav_file, index_col=0)

    for ch in channels:
        if ch not in fiberbehav_df.columns:
            continue
        for cfg_name, cfg in configs.items():
            record, events_df = quantif.extract_transients(
                fiberbehav_df, mouse, batch, group,
                zone_cols=None, behav_cols=['Sniff Novel', 'Sniff HC', 'Sniff Clean'],
                baseline_col=None,
                baseline_exclude_cols=None,
                merge_into=None,
                stim_sniff_pairs=None,
                dff_col=ch, fps=10, min_distance_s=0.5,
                plot=True,
                plot_path=plot_dir / f"{batch}_{mouse}_{ch.replace(' ', '_')}_{cfg_name}.png",
                **cfg)
            results[(ch, cfg_name)]['summary'].append(record)
            results[(ch, cfg_name)]['events'].append(events_df)

for (ch, cfg_name), res in results.items():
    if not res['summary']:
        continue
    fname = f"transients_{cfg_name}_{ch.replace(' ', '_')}.xlsx"
    with pd.ExcelWriter(repo_path / fname) as writer:
        pd.DataFrame(res['summary']).to_excel(writer, sheet_name='Summary', index=False)
        pd.concat(res['events'], ignore_index=True).to_excel(writer, sheet_name='Events', index=False)

print(f"\n✅ Transient analysis complete. Data saved in: {repo_path}")

# %%
