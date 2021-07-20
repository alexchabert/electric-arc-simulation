# -*- coding: utf-8 -*-
"""
Created on Sun May  2 15:47:47 2021

@author: jonathan.andrea & alexis.chabert
"""

#Reset environment variables
#--------------------------------------------------------------
# from IPython import get_ipython
# get_ipython().magic('reset -sf')

#Libraries
#--------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as inter

#Static characteristic function
#--------------------------------------------------------------
def F(It,a,b,Rc):  
    return a*Rc*It/(Rc*It*np.arctan(b*It)+a)

#Fonction Tau(I) 
#to simulate the time constant during non-arc discharge regimes
#--------------------------------------------------------------
def Tau(I,Tau_min,Tau_max,Seuil_I,Slope):
    return Tau_max*(1/np.pi*np.arctan(Slope*(np.abs(I)-Seuil_I))+0.5)+Tau_min

#Main program
#--------------------------------------------------------------
    
#Simulation parameters
scale=100                           #used to keep a constant number of AC periods
N=25000*scale                       #Length of time vector.
Te=1e-5/scale                       #Sampling period [s].
t=Te*np.linspace(0, N-1, N)         #Time vector

#Circuit parameters
R=13.6                              #Resistor value

#Load experimental files
#--------------------------------------------------------------
path=''
fileIarc=np.loadtxt(path+'Resistance_16A/I.txt')     
fileVg=np.loadtxt(path+'Resistance_16A/V.txt')
fileVarc=np.loadtxt(path+'Resistance_16A/Va.txt')

Iarcexp=fileIarc[:,1]
Varcexp=fileVarc[:,1]
texp=fileVarc[:,0]

Vg_raw=fileVg[:,1][:-2]
t_raw=fileVg[:,0][:-2]
fVg=inter.interp1d(t_raw,Vg_raw,kind='cubic',axis=-1)
Vg=fVg(t)

Iexp=Iarcexp[:int(len(texp)/4)]
Vexp=Varcexp[:int(len(texp)/4)]
tex=0.10001+texp[:int(len(texp)/4)]

#Arc parameter scenario definition
#--------------------------------------------------------------
start=int(9667.6*scale)
stop=int(20906*scale)

a=np.zeros(N)                       #alpha
b=np.zeros(N)                       #beta
Rc=np.zeros(N)                      #Rc
a[start:]=47.12
b[:stop]=0.8*np.ones(stop)
b[stop:]=0.1*np.ones(N-stop)
Rc[:stop]=1900.*np.ones(stop)
Rc[stop:]=3800.*np.ones(N-stop)

Vg[stop:]=0

#Vector initialization
I=np.zeros(N)                       #Arc current vector
V=np.zeros(N)                       #Arc voltage vector

#Main loop
for k in range(1, N-1):
    m=Te/Tau(I[k-1],0.08e-5,4e-5,0.5,1)
    I[k]=(-V[k-1]+Vg[k-1])/R
    V[k]=1./(1+m)*(m*F(I[k],a[k],b[k],Rc[k])+V[k-1])

#Simulation code ENDS here.    
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------


#OPTIONAL:
#Three parts:
#   -Error signal power calculation
#   -Voltage plot
#   -Current plot    
#   -Plot of the evolution of the arc model parameters' 
    
#Selecting arc fault data (croping unnecessary samples)
#--------------------------------------------------------------
Iexp=Iexp[int(len(tex)/4):]
Vexp=Vexp[int(len(tex)/4):]
tex=tex[:int(len(tex)*3/4)]
I=I[int(len(t)/4):]
V=V[int(len(t)/4):]
t=t[:int(len(t)*3/4)]

#Simulation results
#--------------------------------------------------------------
#--------------------------------------------------------------
#Error Power calculation
#--------------------------------------------------------------
errPV=1./len(tex) * np.sqrt(np.sum((V[::10]-Vexp[:])**2))
errPI=1./len(tex) * np.sqrt(np.sum((I[::10]-Iexp[:])**2))
print('error power on voltage =',errPV)
print('error power on current =',errPI)

#Plots of the simulation results 
#--------------------------------------------------------------
#Parameters of the plots
#--------------------------------------------------------------
plt.close('all')
plt.rc('font', family='DejaVu Sans', serif='Times')
plt.rc('text', usetex=False)
plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=10)
plt.rc('axes', labelsize=10)

#Plot voltages
#--------------------------------------------------------------
fig, axs = plt.subplots(2, 1,sharex='col')
axs[0].set_xlim([0,t[-1]*1e3])

axs[0].plot(t*1e3, V,'tab:blue')
axs[0].set(xlabel=r'Time (ms)',ylabel=r'Sim. arc voltage (V)')
axs[0].get_yaxis().set_label_coords(-0.15,0.5)

axs[0].annotate('', xy=((start-N*1/4)*Te*1e3,-200), 
                xytext=((stop-N*1/4)*Te*1e3,-200), 
                arrowprops=dict(arrowstyle='<->'))
axs[0].annotate('Arc',
            xy=(1, 0), xycoords='axes fraction',
            xytext=(0.5, 0.05), textcoords='axes fraction',
            horizontalalignment='right',
            verticalalignment='bottom')

axs[1].plot(tex*1e3, Vexp, 'tab:blue')
axs[1].set(xlabel=r'Time (ms)',ylabel=r'Exp. arc voltage (V)')
axs[1].get_yaxis().set_label_coords(-0.15,0.5)

axs[0].set_ylim([-220,220])
axs[1].set_ylim([-220,220])

fig.subplots_adjust(wspace=0.01,hspace=0.1)
width  = 3.487
height = width / 1.618 *2
fig.set_size_inches(width, height)

for ax in fig.get_axes():
    ax.label_outer()
    
axs[0].axvline(x=34,color='k',linewidth=0.25,linestyle='--')     
axs[1].axvline(x=34,color='k',linewidth=0.25,linestyle='--')  

axs[0].axvline(x=147,color='k',linewidth=0.25,linestyle='--')     
axs[1].axvline(x=147,color='k',linewidth=0.25,linestyle='--')  

#Plot current
#--------------------------------------------------------------
fig2, axs = plt.subplots(2, 1,sharex='col')
axs[0].set_xlim([0,t[-1]*1e3])

axs[0].plot(t*1e3, I, 'tab:green')
axs[0].set(ylabel=r'Sim. arc current (A)')
axs[0].get_yaxis().set_label_coords(-0.15,0.5)

axs[0].annotate('', xy=((start-N*1/4)*Te*1e3,-28), 
                xytext=((stop-N*1/4)*Te*1e3,-28), 
                arrowprops=dict(arrowstyle='<->'))
axs[0].annotate('Arc',
            xy=(1, 0), xycoords='axes fraction',
            xytext=(0.5, 0.05), textcoords='axes fraction',
            horizontalalignment='right',
            verticalalignment='bottom')

axs[1].plot(tex*1e3, Iexp, 'tab:green')
axs[1].set(xlabel=r'Time (ms)',ylabel=r'Exp. arc current (A)')
axs[1].get_yaxis().set_label_coords(-0.15,0.5)

axs[0].set_ylim([-30,30])
axs[1].set_ylim([-30,30])

fig2.subplots_adjust(wspace=0.01,hspace=0.1)
width  = 3.487
height = width / 1.618 *2
fig2.set_size_inches(width, height)

for ax in fig2.get_axes():
    ax.label_outer()
      
axs[0].axvline(x=34,color='k',linewidth=0.25,linestyle='--')     
axs[1].axvline(x=34,color='k',linewidth=0.25,linestyle='--')  
axs[0].axvline(x=147,color='k',linewidth=0.25,linestyle='--')     
axs[1].axvline(x=147,color='k',linewidth=0.25,linestyle='--')  
  
#Arc model parameters plot
#--------------------------------------------------------------
#Selecting the relevant data
#--------------------------------------------------------------
a=a[int(N/4):]
b=b[int(N/4):]
Rc=Rc[int(N/4):]
#Plots of the evolution of the arc model parameters'
#--------------------------------------------------------------
fig3, (ax1, ax2,ax3) = plt.subplots(3, sharex=True)
ax1.set_xlim([0,t[-1]*1e3])

ax1.plot(t*1e3, a)
ax1.set_ylabel(r'$\alpha$''\n (V)',rotation=0)
ax1.get_yaxis().set_label_coords(-0.2,0.25)
ax1.set_ylim([-10,100])
ax1.axvline(x=34,color='k',linewidth=0.25,linestyle='--')  
ax1.axvline(x=147,color='k',linewidth=0.25,linestyle='--') 

ax1.annotate('', xy=((start-N*1/4)*Te*1e3,70), 
                xytext=((stop-N*1/4)*Te*1e3,70), 
                arrowprops=dict(arrowstyle='<->'))
ax1.annotate('Arc',
            xy=(1, 0), xycoords='axes fraction',
            xytext=(0.5, 0.75), textcoords='axes fraction',
            horizontalalignment='right',
            verticalalignment='bottom')

ax1.annotate('', xy=(Te*1e3,70),
             xytext=((start-N*1/4)*Te*1e3,70),
             arrowprops=dict(arrowstyle='<->'))
ax1.annotate('CC',
            xy=(0.09, 0.75), xycoords='axes fraction',
            xytext=(0.09, 0.75), textcoords='axes fraction',
            horizontalalignment='center',
            verticalalignment='bottom')

ax1.annotate('', xy=((stop-N*1/4)*Te*1e3,70),
             xytext=(3/4*N*Te*1e3,70),
             arrowprops=dict(arrowstyle='<->'))
ax1.annotate('OC',
            xy=(0.89, 0.75), xycoords='axes fraction',
            xytext=(0.89, 0.75), textcoords='axes fraction',
            horizontalalignment='center',
            verticalalignment='bottom')

ax2.plot(t*1e3, b)
ax2.get_yaxis().set_label_coords(-0.2,0.25)
ax2.set_ylabel(r'$\beta$'' \n (a.u.)',rotation=0)
ax2.set_ylim([0,1])
ax2.axvline(x=34,color='k',linewidth=0.25,linestyle='--')  
ax2.axvline(x=147,color='k',linewidth=0.25,linestyle='--') 

ax3.plot(t*1e3, Rc)
ax3.get_yaxis().set_label_coords(-0.2,0.25)
ax3.set_ylabel('Rc\n ($\Omega$)',rotation=0)
ax3.set_ylim([1800,4000])
ax3.set_xlabel(r'Time (ms)')
ax3.axvline(x=34,color='k',linewidth=0.25,linestyle='--')  
ax3.axvline(x=147,color='k',linewidth=0.25,linestyle='--') 

fig3.subplots_adjust(hspace=0.3)
width  = 3.487
height = width / 1.618 *1.5
fig3.set_size_inches(width, height)

plt.show()
