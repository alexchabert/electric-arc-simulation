# -*- coding: utf-8 -*-
"""
Created on Sun May  2 11:46:36 2021

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

#Static characteristic function
#--------------------------------------------------------------
def F(It):  
    a=47
    Rc=3000
    b=1.47  
    return a*Rc*It/(Rc*It*np.arctan(b*It)+a)

#Main program
#--------------------------------------------------------------
#Simulation parameters
#--------------------------------------------------------------
N=1000000                            #Number of points in the simulation
Te=1e-7                             #Sampling period
t=Te*np.linspace(0, N-1, N)         #Time vector

#Circuit parameters
#--------------------------------------------------------------
R=14                                #Resistor value
L=3e-3                              #Inductor value
Tau=1e-5                            #Arc time constant
f=50                                #Generator frequency
Vgmax=300                           #Generator amplitude
Vg=Vgmax*np.sin(2*np.pi*f*t)        #Generator signal

#ABCD matrix parameters
#--------------------------------------------------------------
m=Te/Tau                            #Constant to simplify the code
A= np.array([[1-Te*R/L, -Te/L], 
             [0, 1/(1+m)]])         #A matrix
B= np.array([[Te/L, 0], 
             [0, m/(1+m)]])         #B matrix

#Initialisation of the state space vectors
#--------------------------------------------------------------
X=np.zeros((2,N))                   #State space Vector X
U=np.zeros((2,1))                   #Input vector

#Main loop
#--------------------------------------------------------------
for k in range(1, N-1):
    U[0,0]=Vg[k-1]
    U[1,0]=F(X[0,k-1]+(1/m)*(X[0,k-1]-X[0,k-2]))
    X[:,k]=np.dot(A,X[:,k-1])+np.dot(B,U).transpose()

#Raw plot simulation results
#--------------------------------------------------------------
plt.close()
plt.plot(t,Vg)                      #Generator voltage
plt.plot(t,np.array(X[0,:]))        #Current
plt.plot(t,np.array(X[1,:]))        #Voltage