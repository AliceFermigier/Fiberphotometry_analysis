# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 11:41:03 2021

Plot tracking

@author: Alice Fermigier
"""

import matplotlib.pyplot as plt
import pandas as pd

df_tracking = pd.read_csv('E:/Alice/Fiber/202103_CDHFDdvHPC/20210302_ORM/Test1h/CD1_tracking.csv')

fig3 = plt.figure(figsize=(20,20))
ax2 = fig3.add_subplot(111)

p1, = ax2.plot('Item2.X', 'Item2.Y', linewidth=.5, color='black', data=df_tracking)