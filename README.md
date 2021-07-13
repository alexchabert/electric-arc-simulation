# electric-arc-simulation
Method to simulate Andrea's electric arc model with a state-space formalism described in [cite publication][submitted to IEEE].

This repository contains three arc simulation codes described as follows.

--------------------------------------------Data Description-----------------------------------------------------------

'R_Arc_AC_git': Arc simulation made to fit an experimental arc in series with an R load under an AC generator.

'RL_Arc_DC_git': Arc simulation made to fit an experimental arc in series with an RL load under a DC generator. 
Ignition is done by contact opening.

'RL_Arc_ABCD_git': Simulation of an arc in series with an RL load under an AC generator. 
The purpose of this code is to show the explicit definition of the ABCD matrices.

---------------------------------------------General Infos-------------------------------------------------------------

'R_AC' and 'RL_DC' folders allow to reproduce the curves presented in [cite publication][submitted to IEEE].

Loading paths for experimental data are relative (data should be kept in the code repertory).

--------------------------------------------Libraries & Test-----------------------------------------------------------

The libraries required are: 

  -numpy
  
  -matplotlib
  
  -scipy (to interpolate generator data)
  
  -pandas (to open the .csv file)

'environment.yml' contains the Anaconda environment with these libraries.

The codes have been tested on Python 3.7. 


  
