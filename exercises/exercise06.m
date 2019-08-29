%% Solution by Stefan Sicklinger 2019: CA-CFAR
clc;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data_points
Ns = 1000;
% Generate random noise
s=abs(randn(Ns,1));
%Targets location. Assigning bin 100, 200, 300 and 700 as Targets
%with the amplitudes of 8, 9, 4, 11.
s([100 ,200, 300, 700])=[8 9 4 11];
%plot the output
plot(s);
% TODO: Apply CFAR to detect the targets by filtering the noise.
%Training Cells
T=12;
%Guard Cells
G=4;
% Offset : Adding room above noise threshold for desired SNR
offset=5;
% Vector to hold threshold values
threshold_cfar = [];
%Vector to hold final signal after thresholding
signal_cfar = [];
% 2. Slide window across the signal length
for i = 1:(Ns-(G+T+1))
    noise_level=sum(s(i:i+T-1));
    threshold= (noise_level/T)*offset;
    threshold_cfar =[threshold_cfar,{threshold}];
    signal =s(i+T+G);
    % 8. Filter the signal above the threshold
    if(signal<threshold)
        signal=0;
    end
    signal_cfar = [signal_cfar, {signal}];
end
% plot the filtered signal
plot (cell2mat(signal_cfar),'g--');
% plot original sig, threshold and filtered signal within the same figure.
figure,plot(s);
hold on,plot(cell2mat(circshift(threshold_cfar,G)),'r--','LineWidth',2)
hold on, plot (cell2mat(circshift(signal_cfar,(T+G))),'g--','LineWidth',4);
legend('Signal','CFAR Threshold','detection')