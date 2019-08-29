%% Solution by Stefan Sicklinger 2019
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Operating frequency (Hz)
fc = 77.0e9;
%Transmitted power (W)
Pt = 3e-3;
%Antenna Gain (linear) (dBi)
G =  10000;
%Minimum Detectable Power (dBm)
Ps = 1e-10;
%RCS of a car (m^2)
RCS = 100;
%Speed of light (m/s)
c = 3*10^8;
%Calculate the wavelength (m)
lambda = c/fc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Measure the maximum Range a Radar can see. (m)
numerator=Pt*G*G*lambda*lambda*RCS;
denominator=Ps*(4*pi)^3;
R=nthroot(numerator/denominator,4);
disp(R);