# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 11:24:30 2022

Code to vizualize data in Dash

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
import Fiberpho_process

#%%

exp = 'OdDis post shock'
session = 'Test'
mouse = 'B1m'
data_type = 'fiberpho'

exp_path = analysis_path / f'{exp}'
session_path = exp_path / f'{session}'
data_path_exp = data_path / '20230222_AliceF_CA2b5OdDispostshock'

rawdata_path = data_path_exp / f'{mouse}_1.csv'

app = Dash(__name__)

if data_type == 'fiberpho':
    deinterleaved_df = Fiberpho_process.deinterleave(pd.read_csv(rawdata_path, skiprows=1, usecols=['Time(s)','AIn-1','DI/O-1','DI/O-2']))
    fig = px.line(deinterleaved_df[50:], x='Time(s)', y='470 Deinterleaved')
    #fig = px.line(deinterleaved_df, x='Time(s)', y='405 Deinterleaved')
    
elif data_type == 'plethysmo':
    plethys_df = pd.read_csv(rawdata_path, skiprows=1, usecols=['Time(s)','AIn-4'])
    plethys_df = plethys_df.loc[[i for i in range(0,len(plethys_df),600)]] #downsample plethys_df
    fig = px.line(plethys_df, x='Time(s)', y='405 Deinterleaved')

app.layout = html.Div([
    html.H4(f'{exp} {session} {mouse}'),
    dcc.Graph(
        id=f'{exp}',
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

# fiberpho = pd.read_csv(data_path_exp / f'{mouse}_1_dFFfilt.csv')

# plt.plot('Time(s)', 
#           'Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0-Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass', 
#           data = fiberpho1)

# plt.plot('Time(s)', 
#           'Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0-Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass', 
#           data = fiberpho2)

# plt.plot('Time(s)', 
#           'Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0-Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass', 
#           data = fiberpho)

# plt.plot('Time(s)', 
#           'AIn-4', 
#           data = rawdata)

# fiberpho.to_csv(data_path_exp / f'{mouse}_1_dFFfilt.csv', index=False)
# #rawdata.to_csv(data_path_exp / f'{mouse}_1.csv')

# fiberpho.drop(columns='Unnamed: 0.1', inplace = True)


#for interactive plotting, maybe later
"""
def visualize(filepath, x, y1, y2, y3, plot_name, removeArtifacts):
    

	# plotting control and signal data

	if (y1==0).all()==True:
		y1 = np.zeros(x.shape[0])

	coords_path = os.path.join(filepath, 'coordsForPreProcessing_'+plot_name[0].split('_')[-1]+'.npy')
	name = os.path.basename(filepath)
	fig = plt.figure()
	ax1 = fig.add_subplot(311)
	line1, = ax1.plot(x,y1)
	ax1.set_title(plot_name[0])
	ax2 = fig.add_subplot(312)
	line2, = ax2.plot(x,y2)
	ax2.set_title(plot_name[1])
	ax3 = fig.add_subplot(313)
	line3, = ax3.plot(x,y2)
	line3, = ax3.plot(x,y3)
	ax3.set_title(plot_name[2])
	fig.suptitle(name)

	hfont = {'fontname':'Helvetica'}

	ax3.set_xlabel('Time(s)', **hfont)

	global coords
	coords = []

	# clicking 'space' key on keyboard will draw a line on the plot so that user can see what chunks are selected
	# and clicking 'd' key on keyboard will deselect the selected point
	def onclick(event):
		#global ix, iy

		if event.key == ' ':
			ix, iy = event.xdata, event.ydata
			print('x = %d, y = %d'%(
			    ix, iy))

			y1_max, y1_min = np.amax(y1), np.amin(y1)
			y2_max, y2_min = np.amax(y2), np.amin(y2)

			#ax1.plot([ix,ix], [y1_max, y1_min], 'k--')
			#ax2.plot([ix,ix], [y2_max, y2_min], 'k--')

			ax1.axvline(ix, c='black', ls='--')
			ax2.axvline(ix, c='black', ls='--')
			ax3.axvline(ix, c='black', ls='--')

			fig.canvas.draw()

			global coords
			coords.append((ix, iy))

			#if len(coords) == 2:
			#    fig.canvas.mpl_disconnect(cid)

			return coords

		elif event.key == 'd':
			if len(coords)>0:
				print('x = %d, y = %d; deleted'%(
			    	coords[-1][0], coords[-1][1]))
				del coords[-1]
				ax1.lines[-1].remove()
				ax2.lines[-1].remove()
				ax3.lines[-1].remove()
				fig.canvas.draw()

			return coords

	# close the plot will save coordinates for all the selected chunks in the data
	def plt_close_event(event):
		global coords
		if coords and len(coords)>0:
			name_1 = plot_name[0].split('_')[-1]
			np.save(os.path.join(filepath, 'coordsForPreProcessing_'+name_1+'.npy'), coords)
			print('Coordinates file saved at {}'.format(os.path.join(filepath, 'coordsForPreProcessing_'+name_1+'.npy')))
		fig.canvas.mpl_disconnect(cid)
		coords = []


	cid = fig.canvas.mpl_connect('key_press_event', onclick)
	cid = fig.canvas.mpl_connect('close_event', plt_close_event)
	#multi = MultiCursor(fig.canvas, (ax1, ax2), color='g', lw=1, horizOn=False, vertOn=True)

	#plt.show()
	#return fig

# function to plot control and signal, also provide a feature to select chunks for artifacts removal
def visualizeControlAndSignal(filepath, removeArtifacts):
	path_1 = find_files(filepath, 'control_*', ignore_case=True) #glob.glob(os.path.join(filepath, 'control*'))
	
	path_2 = find_files(filepath, 'signal_*', ignore_case=True) #glob.glob(os.path.join(filepath, 'signal*'))
	

	path = sorted(path_1 + path_2, key=str.casefold)
	
	if len(path)%2 != 0:
		raise Exception('There are not equal number of Control and Signal data')
	
	path = np.asarray(path).reshape(2,-1)
	
	for i in range(path.shape[1]):
		
		name_1 = ((os.path.basename(path[0,i])).split('.')[0]).split('_')
		name_2 = ((os.path.basename(path[1,i])).split('.')[0]).split('_')
		
		ts_path = os.path.join(filepath, 'timeCorrection_'+name_1[-1]+'.hdf5')
		cntrl_sig_fit_path = os.path.join(filepath, 'cntrl_sig_fit_'+name_1[-1]+'.hdf5')
		ts = read_hdf5('', ts_path, 'timestampNew')
		
		control = read_hdf5('', path[0,i], 'data').reshape(-1)
		signal = read_hdf5('', path[1,i], 'data').reshape(-1)
		cntrl_sig_fit = read_hdf5('', cntrl_sig_fit_path, 'data').reshape(-1)

		plot_name = [(os.path.basename(path[0,i])).split('.')[0], 
					 (os.path.basename(path[1,i])).split('.')[0],
					 (os.path.basename(cntrl_sig_fit_path)).split('.')[0]]
		visualize(filepath, ts, control, signal, cntrl_sig_fit, plot_name, removeArtifacts)

"""


