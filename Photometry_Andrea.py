#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# In[2]:


os.getcwd()


# In[3]:


os.chdir('C:\\Users\\Andrea\\Desktop\\fiber_photometry')


# In[4]:


data = pd.read_csv ('79157_pfc_basal.csv')


# In[5]:


# Preprocessing of data. I take the data starting at 150 s

data.columns = data.loc[0]
data = data.drop(0)
data.columns = ['Time', 'ISOS', 'Dlight', 'Raw', 'Out1', 'Out2', 'NaN']
for col in data.columns: 
    data[col] = pd.to_numeric(data[col])
data = data.drop('NaN', axis = 1)
data = data.dropna()
data = data[data['Time'] > 150.0]
data


# In[45]:


# extractin column data into different array

ISOS = np.array(data['ISOS'])
Dlight = np.array(data['Dlight'])
raw = np.array(data['Raw'])
time = np.array(data['Time'])

# I found on doric studio that the acquisition rate was set at 12Ksps
sampling_rate = 12000


# In[46]:


# plot of raw data to see what the signal looks like

plt.figure(figsize=(10, 8))
plt.plot (time, ISOS, label ='isos')
plt.plot (time, Dlight, 'r', label = 'dlight')
plt.legend()


# In[12]:


# downloading scipy modules

from scipy.signal import medfilt, butter, filtfilt
from scipy.stats import linregress


# In[13]:


# I found on the internet that people use a low pass filter to reduce noise signal

b, a = butter(2, 5, btype='low', fs = sampling_rate)
dlight_denoised = filtfilt(b,a, Dlight)
isos_denoised = filtfilt(b, a, ISOS)


# In[14]:


plt.figure(figsize=(12, 12))
plt.plot(time, dlight_denoised, label = 'dlight_denoised')
plt.plot(time, isos_denoised, label ='ISOS_denoised')
plt.title ('Denoised signals')

plt.legend()


# In[15]:


# Application of a highpass filter for photobleaching correction. Found on the internet... 

b,a = butter(2, 0.001, btype='high', fs =sampling_rate)
dlight_highpass = filtfilt(b, a, dlight_denoised)
isos_highpass = filtfilt(b, a, isos_denoised)


# In[16]:


plt.figure(figsize=(12, 12))
plt.plot(time, dlight_highpass,'g', label = 'dlight highpass')
plt.plot(time, isos_highpass - 0.005, 'r', label = 'isos highpass')
plt.xlabel ('Time(s)')
plt.ylabel ('Signal (Volts)')
plt.title ('Photobleaching correction')
plt.legend()


# In[19]:


# In github I found a code for calculating the dF / F based on motion correction of signals after performing a regression analysis of the isos and dlight signal
# I do not understand comppletely this concept, but I tried to apply it to my data

slope, intercept, r_value, p_value, std_err = linregress(x=isos_highpass, y= dlight_highpass)
plt.scatter (isos_highpass[::5], dlight_highpass[::5], alpha= 0.1, marker ='.')
x = np.array(plt.xlim())
plt.plot(x, intercept+slope*x)


# In[22]:


dlight_est_motion = intercept + slope * isos_highpass
dlight_corrected = dlight_highpass - dlight_est_motion
plt.plot(time, dlight_corrected, label = 'corrected')
plt.plot(time, dlight_highpass,'g', label = 'not-corrected')
plt.plot(time, dlight_est_motion, 'y', label = 'estimated motion')


# In[23]:


b, a = butter(2, 0.001, btype = 'low', fs = sampling_rate)
baseline_fluorescence = filtfilt(b, a, dlight_denoised, padtype='even')
plt.figure(figsize=(12,12))
plt.plot(time, dlight_denoised, 'g', label = 'dlight_denoised')
plt.plot(time, baseline_fluorescence, 'k', label = 'baseline fluorescence ')


# In[24]:


dlight_dF_F = dlight_corrected / baseline_fluorescence
plt.figure(figsize=(20, 8))
plt.plot(time, dlight_dF_F)
plt.xlabel ('Time(s)')
plt.ylabel ('dlight dF-F (%)')


# In[25]:


# plotting the graph on amphetamine 20 min 

data2 = pd.read_csv ('79157_pfc_amp_20 min.csv')
data2.columns = ['Time', 'ISOS', 'Dlight', 'Raw', 'Out1', 'Out2', 'NaN']
data2 = data2.drop('NaN', axis = 1)
# removing first row
data2 = data2.drop(0)


# In[26]:


#converting columns values in numbers 
for col in data2.columns:
    data2[col] = pd.to_numeric(data2[col])


# In[27]:


# removing NaN
data2 = data2.dropna()

# removing all the time points lower than 400 seconds to avoid electrical artifact. 
# Do you know any filter to remove electrical artifact ? 
# I tried to use the scipy median filter but it doesn't work. Any suggestions ? 

data2 = data2[data2['Time'] > 400]
data2


# In[28]:


dlight2 = np.array(data2['Dlight'])
isos2 = np.array(data2['ISOS'])
time2 = np.array(data2['Time'])
sampling_rate = 12000


# In[29]:


plt.figure(figsize=(12,12))
plt.plot(time2, dlight2, label = 'dlight')
plt.plot (time2, isos2, label = 'isos')
plt.legend()


# In[30]:


#applying  a filter to reduce noise 
b, a = butter(2, 5, btype = 'low', fs = sampling_rate)
dlight2_denoised = filtfilt(b, a, dlight2)
isos2_denoised = filtfilt(b, a, isos2)


# In[33]:


# Plotting raw vs denoised signal
plt.figure(figsize=(20,8))
plt.title ('Effect of low pass filter')
plt.plot (time2, dlight2, label = 'dlight raw')
plt.plot (time2, dlight2_denoised,label =' dlight_denoised')
plt. legend()


# In[37]:


#photobleaching correction 
b, a = butter(2, 0.001, btype = 'high', fs = sampling_rate)
dlight2_highpass = filtfilt(b, a, dlight2_denoised)
isos2_highpass = filtfilt (b, a, isos2_denoised)


# In[38]:


# Motion correctionfor determining dF/F

slope, intercept, r_value, p_value, sts_err = linregress (x = isos2_highpass, y=dlight2_highpass)
dlight2_est_motion = intercept + slope * isos2_highpass
dlight2_corrected = dlight2_highpass - dlight2_est_motion


# In[39]:


b, a = butter(2, 0.001, btype= 'low', fs = sampling_rate)
baseline_fluorescence = filtfilt(b, a, dlight2_denoised)


# In[40]:


dlight2_dF_F = dlight2_corrected/baseline_fluorescence 


# In[42]:


plt.figure(figsize=(20, 8))
plt.title ('dlight dF_F')
plt.plot (time2, dlight2_dF_F * 100, label = 'dlight2')
plt.xlabel ('Time(s)')
plt.ylabel ('dlight dF_F (%)')
plt.legend()


# In[44]:


# Plotting only 20 second of signal for both basal and stimulated condition

fig, axs = plt.subplots (2, figsize=(12, 12))
fig.tight_layout(pad=3.0)
axs[0].plot(time, dlight_dF_F, label = ' basal condition')
axs[1].plot (time2, dlight2_dF_F, label = 'amphetamine 5 mg/kg')
axs[0].set_title('Basal Condition')
axs[1].set_title ('amphetamine 5 mg/kg')
axs[0].set_xlim(280, 300)
axs[1].set_xlim(400, 420)


# In[ ]:




