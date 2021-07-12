# electric-arc-simulation
New method to simulate Andrea's electric arc model with a state-space formalism described in [cite publication].

This repository contains three arc simulation codes described as follows.

'R_Arc_AC_git': Arc simulation made to fit an experimental arc in series with an R load under an AC generator

'RL_Arc_DC_git': Arc simulation made to fit an experimental arc in series with an RL load under a DC generator. Ignition is done by contact opening

'RL_Arc_ABCD_git': Simulation of an arc in series with an RL load under an AC generator. The purpose of this code is to show the explicit definition of the ABCD matrices.

'R_AC' and 'RL_DC' folders allow to reproduce the curves presented in [cite publication].

Loading paths for experimental data are relative (data should be kept in the code repertory).

The libraries required are: 
  -numpy
  -matplotlib
  -scipy (to interpolate generator data)
  -pandas (to open the .csv file)

The codes have been tested on Python 3.7. 
  
