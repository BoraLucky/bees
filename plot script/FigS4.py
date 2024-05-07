# -*- coding: utf-8 -*-

import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import numpy as np

import matplotlib.pyplot as plt
#from brokenaxes import brokenaxes
import numpy as np
labels = ['0-50bp','50-100bp', '100-1000bp', '1000-10000bp', '10000+bp']  
Del = [11833,2434,557,87,6]  
Dup = [3414,138,111,5,1]     
Inv = [2299,62,55,7,1]      
x = np.arange(len(labels))  
width = 0.28 

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
fig.subplots_adjust(hspace=0.05)  # adjust space between axes


rects1 = ax2.bar(x , Del, width, label='Del')  
rects2 = ax2.bar(x + width, Dup, width, label='Dup')
rects3 = ax2.bar(x + 2*width, Inv, width, label='Inv')
ax2.set_ylabel('Numbers')  
ax2.set_xticks(x+width)      
ax2.set_xticklabels(labels)  
def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax2.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

autolabel(rects1)
autolabel(rects2)
autolabel(rects3)


rects1 = ax1.bar(x , Del, width, label='Del') 
rects2 = ax1.bar(x + width, Dup, width, label='Dup')
rects3 = ax1.bar(x + 2*width, Inv, width, label='Inv')
ax1.set_title('Numbers by lenth and SV_type') 
ax1.legend() 
def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax1.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
autolabel(rects1)
autolabel(rects2)
autolabel(rects3)


ax1.set_ylim(400, 14000) 
ax1.spines['bottom'].set_visible(False)


ax2.set_ylim(0, 150) 
ax2.spines['top'].set_visible(False)

ax1.xaxis.tick_top()
ax1.tick_params(labeltop=False) 
ax2.xaxis.tick_bottom()


d = .5  
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)

fig.tight_layout()
plt.savefig('SV_bar.pdf',dpi=300)
