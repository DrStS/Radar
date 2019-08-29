%% Solution by Stefan Sicklinger 2019
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% beat frequencies (MHz)
fb=[0 1.1 13 24]*1e6;
%The radar maximum range = 300m
Rmax= 300;
%The range resolution = 1m
dres=1;
%Speed of light (m/s)
c = 3*10^8;
%Calculate the wavelength (m) 
lambda = c/fc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the Bsweep of chirp for 1 m resolution (Hz)
Bs=c/(2*dres);
% Calculate the chirp time based on the Radar's Max Range
Ts=5.5*2*Rmax/c;
%define the frequency shifts 
calculated_range=c*Ts*fb/(2*Bs);
% Display the calculated range
disp(calculated_range);