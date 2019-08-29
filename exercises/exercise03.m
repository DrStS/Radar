%% Solution by Stefan Sicklinger 2019
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% targets with following doppler frequency shifts (Hz)
% minus means traveling form the ego car -> receding
% plus means traveling to the ego car -> approaching
deltaFdoppler=[3 , -4.5 , 11 , -3 ]*1e3;
% frequencies (Hz)
frequency = 77e9;   
%Speed of light (m/s)
c = 3*10^8;
%Calculate the wavelength (m)
lambda = c/frequency;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
velociy_targets= deltaFdoppler*lambda/2;
disp(velociy_targets);
