echo on
clc

%                        LIST OF PWATOOL COMMANDS
%
%
%    (1) pwacreate(n,m, 'myfilename.m')                        Approximates a nonlinear 
%                                                              model with a PWADIs 
% 
%
%    (2) pwacreate(n,m, 'myfilename.m', NR)                    Imports a PWA model
% 
%
%
%    (3) y=pwaanalysis(pwasys, setting)                        Checks the stability of a 
%        y=pwaanalysis(pwasys, setting, gain)                  PWA or nonlinear model at
%                                                              a given point xcl.
%                                                              'y' is set to 1 if pwasys is stable.
%                                                              'gain' is the controller Kbar
%                                                   
%
%    (4) [Control, DataSim]=pwasynth(pwasys, x0, setting)       Synthesizes PWA controllers
%                                                               for a PWA or nonlinear model
%                                                               'Control' contains the gains
%                                                               'DataSim' contains the simulation data from Simulink

%
pause; %strike any key to continue
echo off
 


