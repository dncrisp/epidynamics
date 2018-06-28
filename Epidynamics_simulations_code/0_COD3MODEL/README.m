% This folder contains scripts to simulate the model described in Saggio et
% al. 2017

% To produce a simulation from the model use 'Integrator_Cod3_Model.m'
% This script uses: the model in 'Cod3_Model.m', 
%                               which uses the functions:
%                               'eval_resting_state_cartesian.m' 
%                               'Parametrization_GreatCircle.m'
%                   the function in 'Parametrization_VectorsGreatCircle.m'

% The great circle produced by different settings for A and B can be tested
% using 'Test_Parametrization_GreatCircle.m'

% 'Transitions_Cod3_Integrator.m' allows to obtain transitions between
% classes. It uses 'Transitions_Cod3_Model.m'

% The folder 'Figures' contains flat maps of the unfolding with information on the
% frequency and amplitude of the stable limit cycle

% The folder 'FOCUS_r0p4' contains spherical maps of the unfolding and
% information on the amplitude and frequency of the stable limit cycle. It
% also contains the bifurcation curves of the unfolding as computed using
% Matcont and CL_Matcont, and a script to plot this unfolding