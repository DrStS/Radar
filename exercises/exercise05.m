%% Solution by Stefan Sicklinger 2019: 2D FFT Radar: Range Doppler Map(RDM)
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = peaks(20);
X = repmat(P,[5 10]);
imagesc(X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO : Compute the 2-D Fourier transform of the data.  
% Shift the zero-frequency component to the center of the output, and 
% plot the resulting 100-by-200 matrix, which is the same size as X.