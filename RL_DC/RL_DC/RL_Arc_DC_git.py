# -*- coding: utf-8 -*-
"""
Created on Sun May  2 18:31:10 2021

@author: jonathan.andrea & alexis.chabert
"""

#Reset environment variables
#--------------------------------------------------------------
#from IPython import get_ipython
#get_ipython().magic('reset -sf')

#Libraries
#--------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#Static characteristic function
#--------------------------------------------------------------
def F(It,a,b,Rc,n):  
    return a*Rc*It/(Rc*It*np.arctan(b*It**n)+a)

#Loading experimental data
#--------------------------------------------------------------
path=''
file=pd.read_csv(path+'arc_moyen_ok.CSV',delimiter=",",skiprows=1)

Iarcexp=file.values[:,1]
texp=file.values[:,0]
Iexp=file.values[:,1]
Vexp=file.values[:,2]
Vgraw=file.values[:,3]

#Main program
#--------------------------------------------------------------

#Simulation parameters
#--------------------------------------------------------------
N=119559                            #Length of time vector.
Te=6e-7                             #Sampling period [s].
t=Te*np.linspace(0, N-1, N)         #Vecteur temps

#Circuit parameters
R=23.42                             #Resistor value
L=1.01e-3                           #Inductor value
Tau=1e-5                            #Arc time constant
v=2200                              #System opening speed

#Use experimental voltage source as simulation voltage source
#--------------------------------------------------------------
Vg=Vgraw

#Arc model parameters' evolution
#--------------------------------------------------------------
start=int(37134)
stop=int(86126)

n=np.ones(len(t))
Len=np.zeros(len(t))
for k in range(len(t)):
    if (k>=start):     
        Len[k]=v*(t[k]-t[start])**2
    n[k]=5.20531228e-01*(1-np.exp(-Len[k]*1.96574187e+01*7.51862417e-03)
              *np.sin(Len[k]*5.36472293e-02+2.53404315e-01)
              /np.sqrt(1-7.51862417e-03**2))
        
a=np.zeros(len(t)) #alpha
b=np.zeros(len(t)) #beta
Rc=np.zeros(len(t)) #Rc

a[start:]=13.29856856 *1.1286
a=a+11.23337645*np.sqrt(Len) /1.1286

b[:stop]=1.63030998/(Len[:stop]+1)+0.24306177
b[stop:]=0.1*np.ones(N-stop)

Rc[:stop]=2221*np.ones(stop)
Rc[stop:]=2221*np.ones(N-stop)
 
#Vector initialization
#--------------------------------------------------------------
I=np.zeros(N)                       #Arc current vector
V=np.zeros(N)                       #Arc voltage vector

I[0]=Vg[0]/R

#Main loop
#--------------------------------------------------------------
m=Te/Tau                            

for k in range(1, N):
    I[k]=Te*(-R/L*I[k-1]-V[k-1]/L+Vg[k-1]/L)+I[k-1]
    V[k]=1/(1+m)*(m*F(I[k]+(1/m)*(I[k]-I[k-1]),a[k],b[k],Rc[k],n[k])+V[k-1])
 
#Simulation code ENDS here.    
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------



#OPTIONAL:
#Three parts:
#   -Error signal power calculation
#   -Voltage and current plot
#   -Plot of the evolution of the arc model parameters' 

#Error Power calculation
#--------------------------------------------------------------
errPV=1./len(texp) * np.sqrt(np.sum((V[::]-Vexp[:])**2))
errPI=1./len(texp) * np.sqrt(np.sum((I[::]-Iexp[:])**2))
print('error power on voltage =',errPV)
print('error power on current =',errPI)

texp=t

#Plots of the simulation results 
#--------------------------------------------------------------
#Parameters of the plots
#--------------------------------------------------------------
plt.close('all')
plt.rc('font', family='serif', serif='Times')
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=10)
plt.rc('axes', labelsize=10)

#Plot voltages and currents
#--------------------------------------------------------------
fig, axs = plt.subplots(2, 1,sharex='col')
axs[0].set_xlim([0,t[-1]*1e3])
axs[0].plot(texp*1e3, Vexp, 'tab:blue',label='Experimental')
axs[0].plot(t*1e3, V,color='k',label='Simulated')
axs[0].set(xlabel=r'Time (ms)',ylabel=r'Arc voltage (V)')
axs[0].get_yaxis().set_label_coords(-0.15,0.5)
axs[0].set_ylim([-5,100])

axs[1].plot(texp*1e3, Iexp, 'tab:green',label='Experimental')
axs[1].plot(t*1e3, I,color='k',label='Simulated')
axs[1].set(xlabel=r'Time (ms)',ylabel=r'Arc current (A)')
axs[1].get_yaxis().set_label_coords(-0.15,0.5)
axs[1].set_ylim([-0.2,3])

axs[1].annotate('', xy=(start*Te*1e3,0.05), xytext=(stop*Te*1e3,0.05), 
                arrowprops=dict(arrowstyle='<->'))
axs[1].annotate('Arc',
            xy=(1, 0), xycoords='axes fraction',
            xytext=(40, 0.25), textcoords='data',
            horizontalalignment='right',
            verticalalignment='bottom')

axs[0].legend()
axs[1].legend()

fig.subplots_adjust(wspace=0.01,hspace=0.1)
width  = 3.487
height = width / 1.618 *2
fig.set_size_inches(width, height)

for ax in fig.get_axes():
    ax.label_outer()

#plot vertical lines
#--------------------------------------------------------------
axs[0].axvline(x=0.0222804*1e3,color='k',linewidth=0.25,linestyle='--')     
axs[1].axvline(x=0.0222804*1e3,color='k',linewidth=0.25,linestyle='--')  

axs[0].axvline(x=0.05167559*1e3,color='k',linewidth=0.25,linestyle='--')     
axs[1].axvline(x=0.05167559*1e3,color='k',linewidth=0.25,linestyle='--')  

#Plot of the evolution of the arc model parameters' 
#--------------------------------------------------------------
fig3, (ax1,ax2,ax3,ax4) = plt.subplots(4, sharex=True)
ax1.set_xlim([0,t[-1]*1e3])

ax1.plot(t*1e3, Len)
ax1.set_ylabel(r'L''\n (mm)',rotation=0)
ax1.get_yaxis().set_label_coords(-0.2,0.25)
ax1.set_ylim([-0.3,6])
ax1.axvline(x=0.0222804*1e3,color='k',linewidth=0.25,linestyle='--')  
ax1.axvline(x=0.05167559*1e3,color='k',linewidth=0.25,linestyle='--') 

ax1.annotate('', xy=(start*Te*1e3,3.8), xytext=(stop*Te*1e3,3.8), 
             arrowprops=dict(arrowstyle='<->'))
ax1.annotate('Arc',
            xy=(0.5, 0.7), xycoords='axes fraction',
            xytext=(0.5, 0.7), textcoords='axes fraction',
            horizontalalignment='center',
            verticalalignment='bottom')

ax1.annotate('', xy=(Te*1e3,3.8), xytext=(start*Te*1e3,3.8), 
             arrowprops=dict(arrowstyle='<->'))
ax1.annotate('CC',
            xy=(0.15, 0.7), xycoords='axes fraction',
            xytext=(0.15, 0.7), textcoords='axes fraction',
            horizontalalignment='center',
            verticalalignment='bottom')

ax1.annotate('', xy=(stop*Te*1e3,3.8), xytext=(len(texp)*Te*1e3,3.8), 
             arrowprops=dict(arrowstyle='<->'))
ax1.annotate('OC',
            xy=(0.86, 0.7), xycoords='axes fraction',
            xytext=(0.86, 0.7), textcoords='axes fraction',
            horizontalalignment='center',
            verticalalignment='bottom')

ax2.plot(t*1e3, a)
ax2.set_ylabel(r'$\alpha$''\n (V)',rotation=0)
ax2.get_yaxis().set_label_coords(-0.2,0.25)
ax2.set_ylim([-5,40])
ax2.axvline(x=0.0222804*1e3,color='k',linewidth=0.25,linestyle='--')  
ax2.axvline(x=0.05167559*1e3,color='k',linewidth=0.25,linestyle='--') 

ax3.plot(t*1e3, b)
ax3.get_yaxis().set_label_coords(-0.2,0.25)
ax3.set_ylabel(r'$\beta$'' \n (a.u.)',rotation=0)
ax3.set_ylim([0,2.1])
ax3.axvline(x=0.0222804*1e3,color='k',linewidth=0.25,linestyle='--')  
ax3.axvline(x=0.05167559*1e3,color='k',linewidth=0.25,linestyle='--') 

ax4.plot(t*1e3, n)
ax4.get_yaxis().set_label_coords(-0.2,0.25)
ax4.set_ylabel(r'n''\n',rotation=0)
ax4.set_ylim([0.38,0.405])
ax4.set_xlabel(r'Time (ms)')
ax4.axvline(x=0.0222804*1e3,color='k',linewidth=0.25,linestyle='--')  
ax4.axvline(x=0.05167559*1e3,color='k',linewidth=0.25,linestyle='--') 

fig3.subplots_adjust(hspace=0.3)
width  = 3.487
height = width / 1.618 *2
fig3.set_size_inches(width, height)

plt.show()
