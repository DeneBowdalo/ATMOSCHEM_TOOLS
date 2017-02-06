import matplotlib.pyplot as plt


fig, axes = plt.subplots(nrows=4, ncols=4,figsize=(15,13))
fig.patch.set_facecolor('white')

c_l = [-20,0,4,9,20]

count = 0
for ax in axes.flat:
    if count < 14:
        pl = ax.scatter([1,2,3,4,5],[1,2,3,4,5],c=c_l,vmin=-6,vmax=6)
    else:
        ax.axis('off')
    count+=1
    
    ax.tick_params(axis='both', which='major', labelsize=18,pad = 7) 
    ax.set_xticks([1,2,3,4,5])
    ax.set_yticks([1,2,3,4,5])
    ax.set_xticklabels(['0.25','0.5','1.0','2.0','4.0'])
    ax.set_yticklabels(['0.25','0.5','1.0','2.0','4.0'])
    
    ax.set_xlim(1,5)
    ax.set_ylim(1,5)
    
plt.tight_layout(pad = 1.5)
fig.subplots_adjust(bottom=0.08)
fig.subplots_adjust(left=0.10)

plt.annotate('ANMVOC',(-2.14,-0.2),xycoords='axes fraction',fontsize=30, va='top')
plt.annotate('NOx',(-4.65,2.6),xycoords='axes fraction',fontsize=30, va='top',rotation=90)

cbar_ax = fig.add_axes([0.58, 0.15, 0.35, 0.06])

#t = [-20,-10,0,10,20]
t = [-6,-3,0,3,6]
t_str = ['-6','-3','0','+3','+6']
cb = fig.colorbar(pl,orientation='horizontal',cax=cbar_ax,ticks=t)
cb.set_ticklabels(t_str)
cb.ax.tick_params(labelsize=24)
cb.set_label('Simulation Rank\n(0 = Best Possible, 96 = Worst Possible)',fontsize=24)

plt.show()
    

