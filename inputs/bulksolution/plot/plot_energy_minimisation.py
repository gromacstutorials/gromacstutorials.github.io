#!/usr/bin/env python
# coding: utf-8

# In[13]:


import copy
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable


# In[14]:


plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino"],
})


# In[15]:


fontsize = 25
font = {'family': 'sans', 'color':  'black', 'weight': 'normal', 'size': fontsize}
myblack = [0.1,0.1,0.1]
myblue = [0, 0.2, 1]
mygrayblue = [0.5, 0.5, 1]
mygray = [0.5, 0.5, 0.5]
myred = [1, 0, 0]
myblue2 = [23/ 255, 63/ 255, 143/ 255]
myshade = list(np.zeros(8))
for i in range(8):
    myshade[i] = [i/32, i/32, i/8]
my_color_1 = np.array([23,63,143])/255
my_color_2 = np.array([60,174,163])/255
my_color_3 = np.array([255,130,85])/255
my_color_4 = np.array([215,0,0])/255
my_color_5 = np.array([150,150,150])/255


# In[16]:


energy = np.loadtxt('epotmin.xvg') # timetep, ke, pe, press


# In[26]:


fig = plt.figure(figsize=(18.5, 6))

ax1 = fig.add_subplot(131)
plt.plot(energy.T[0], energy.T[1], linewidth = 2)

divider = make_axes_locatable(ax1)
ax1.set_xlabel('steps', fontdict=font)
ax1.set_ylabel('potential energy (kJ/mol)', fontdict=font)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.xlim(-50, 1000)
ax1.minorticks_on()
ax1.tick_params('both', length=10, width=2, which='major', direction='in')
ax1.tick_params('both', length=6, width=1.4, which='minor', direction='in')
ax1.xaxis.set_ticks_position('both')
ax1.yaxis.set_ticks_position('both')
ax1.spines["top"].set_linewidth(2)
ax1.spines["bottom"].set_linewidth(2)
ax1.spines["left"].set_linewidth(2)
ax1.spines["right"].set_linewidth(2)
#ax1.set_xticks([0, 0.5, 1, 1.5, 2.0])
#ax1.set_yticks([7, 20, 33, 47, 60])
#labels = ['$1.3$', '$1.4$', '$1.5$', '$1.6$', '$1.7$']
#ax1.set_yticklabels(labels)
minor_locator_y = AutoMinorLocator(2)
ax1.yaxis.set_minor_locator(minor_locator_y)
minor_locator_x = AutoMinorLocator(2)
ax1.xaxis.set_minor_locator(minor_locator_x)
fig.tight_layout()
plt.savefig('epotmin.png', bbox_inches = 'tight', pad_inches = 0.057)  
plt.show()


# In[ ]:




