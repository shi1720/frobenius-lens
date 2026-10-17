import sys,json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'tests'),str(ROOT/'src')]
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from independent import evolve
from frobenius_lens import Experiment,certificate
F=ROOT/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'serif','font.size':10,'axes.labelsize':10,'axes.titlesize':11,'legend.fontsize':9,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
navy='#183d5b'; teal='#078276'; orange='#ba591f'
def save(fig,name):
 fig.savefig(F/(name+'.pdf'),bbox_inches='tight');fig.savefig(F/(name+'.png'),bbox_inches='tight',dpi=180);plt.close(fig)
# Exact analytic landscape, not a fitted statistic.
A=np.zeros((8,128))
for a in range(8):
 for n in range(1,129):
  b=0;u=n
  while u%2==0:b+=1;u//=2
  A[a,n-1]=min(9,n//2**min(a,b))
fig,ax=plt.subplots(figsize=(7,2.9));im=ax.imshow(A,aspect='auto',origin='lower',extent=[.5,128.5,-.5,7.5],cmap='viridis',vmin=0,vmax=9,interpolation='nearest');ax.set(xlabel='Periodic ring length n',ylabel='Snapshot time exponent a (t = 2^a)',xticks=[1,16,32,48,64,80,96,112,128],yticks=range(8));fig.colorbar(im,ax=ax,label='Independent rule constraints (out of 9)',ticks=[0,1,2,4,6,8,9]);ax.scatter([32,64,128],[6,6,6],c='white',edgecolors=orange,s=55,marker='o',linewidths=1.7);ax.text(4,6.6,'radius r = 4; impulse initial state',color='white',fontsize=9);save(fig,'landscape')
# Exact two-rule witness, with direct repeated updates.
r=4;f=[1,0,1,1,0,0,1,0,1];h=[1,0,1,0,0,0,0,0,0];g=[a^b for a,b in zip(f,h)]
fig,axs=plt.subplots(2,2,figsize=(7,4.8),gridspec_kw={'width_ratios':[2.2,1]})
for col,n in enumerate([128,9]):
 for row,rule in enumerate([f,g]):
  hist=np.array([evolve(rule,n,t,(1,),r,2) for t in range(65)])
  ax=axs[row,col];ax.imshow(hist,origin='upper',aspect='auto',interpolation='nearest',cmap=ListedColormap(['#f1f4f5',navy if row==0 else teal]),extent=[-.5,n-.5,64.5,-.5]);ax.axhline(64,color=orange,lw=1.4);ax.set_title(f'{"f" if row==0 else "g=f+h"}, n={n}');ax.set(xlabel='Cell position',ylabel='Time');ax.set_yticks([0,16,32,48,64])
fig.tight_layout();save(fig,'witness')
# Minimal additional geometry.
es=[Experiment(n,64) for n in [32,64,128]]
ns=list(range(1,17));ranks=[certificate(es+[Experiment(n,64)],4)['rank'] for n in ns]
fig,ax=plt.subplots(figsize=(7,2.8));ax.bar(ns,ranks,color=[teal if v==9 else navy for v in ranks],width=.75);ax.axhline(9,color=teal,ls='--',lw=1);ax.axhline(2,color=orange,ls=':',lw=1.2);ax.text(10.5,2.25,'Original design: rank 2',color=orange);ax.annotate('Smallest single-ring repair',xy=(9,9),xytext=(5.8,10.8),arrowprops={'arrowstyle':'->','color':teal},color=teal,fontsize=9);ax.set(xlabel='Additional impulse ring length (snapshot time 64)',ylabel='Joint rule constraints',xticks=ns,ylim=(0,12),yticks=[0,2,4,6,8,9]);save(fig,'repair')
# Seed limitations: exact degrees, not performance estimates.
fig,ax=plt.subplots(figsize=(7,2.8));ns=np.arange(1,65);imp=[]
for n in ns:
 u=int(n)
 while u%2==0:u//=2
 imp.append(min(9,u))
ax.plot(ns,imp,lw=1.2,c=navy,label='Impulse');ax.plot(ns,np.ones(64),lw=2,c=orange,label='All ones');ax.plot(ns,np.zeros(64),lw=2,c=teal,label='All zero');ax.set(xlabel='Periodic ring length n',ylabel='Independent rule constraints',ylim=(-.4,9.8),yticks=[0,1,3,5,7,9]);ax.legend(loc='upper center',ncol=3,frameon=False);save(fig,'seed_effect')
# Pipeline with precise contracts.
fig,ax=plt.subplots(figsize=(7,2));ax.axis('off')
boxes=[(.02,.35,.2,.45,'Declare experiment\n(p, r, n, t, seed)'),(.28,.35,.19,.45,'Pull back\nannihilator Q'),(.53,.35,.17,.45,'Combine\nL = lcm Q'),(.77,.35,.21,.45,'Report rank, fiber,\nwitness, and repair')]
for x,y,w,hh,txt in boxes:
 ax.add_patch(plt.Rectangle((x,y),w,hh,transform=ax.transAxes,facecolor='#edf4f3',edgecolor=teal,lw=1));ax.text(x+w/2,y+hh/2,txt,ha='center',va='center',transform=ax.transAxes,fontsize=9)
for x1,x2 in [(.22,.28),(.47,.53),(.70,.77)]:ax.annotate('',xy=(x2,.575),xytext=(x1,.575),xycoords='axes fraction',arrowprops={'arrowstyle':'->','color':navy})
ax.text(.5,.12,'Independent validation: finite-field row reduction and repeated local updates',ha='center',transform=ax.transAxes,color=navy,fontsize=9);save(fig,'pipeline')
print('Five vector figures generated')
