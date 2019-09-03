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
% initlal image
CFARImage=zeros(1,Ns);
% 1 -> training cell
% 2 -> guard cell
% 3 -> cell under test
CFARImage(1)=1;
imagesc(CFARImage)
colorbar

slideLeftIdx=1;
slideRightIdx=T;
for i = 1:Ns
    noise_level=0;
    lowerIndex=i-G-slideLeftIdx+1;
    if(slideLeftIdx>T)
        slideLeftIdx=T;
    end
    if(lowerIndex>0)
        for k = lowerIndex:(i-G-1)
            CFARImage(k)=1; % coloring training cell
            noise_level=s(k)+noise_level;
        end
        for k = (i-G):(i-1)
            CFARImage(k)=2; % coloring guard cell
        end
        slideLeftIdx=slideLeftIdx+1;
    end
    upperIndex=i+G+slideRightIdx;
    if(slideRightIdx<1)
        slideRightIdx=1;
    end  
    %upperIndex=i+G+T;
    if(upperIndex<=Ns)
        for k = (i+G):(upperIndex)
            CFARImage(k)=1; % coloring training cell
            noise_level=s(k)+noise_level;
        end
        for k = (i+1):(i+G)
            CFARImage(k)=2; % coloring guard cell
        end
    end
    if(upperIndex >= Ns)
       slideRightIdx=slideRightIdx-1;
    end
    %Compute CFAR
    threshold= (noise_level/(2*T))*offset;
    threshold_cfar =[threshold_cfar,{threshold}];
    %Filter the signal above the threshold
    signal =s(i);
    if(s(i)<threshold)
        signal=0;
    end
    signal_cfar = [signal_cfar, {signal}];
    CFARImage(i)=3; % CUT
    imagesc(CFARImage);
    pause(0.0001);
    colorbar;
    CFARImage=zeros(1,Ns);
end
% for i = 1:(Ns-(G+T+1))
%     noise_level=sum(s(i:i+T-1));
%     threshold= (noise_level/T)*offset;
%     threshold_cfar =[threshold_cfar,{threshold}];
%     signal =s(i+T+G);
%     % 8. Filter the signal above the threshold
%     if(signal<threshold)
%         signal=0;
%     end
%     signal_cfar = [signal_cfar, {signal}];
% end
% plot the filtered signal
figure();
plot (cell2mat(signal_cfar),'g--');
% plot original sig, threshold and filtered signal within the same figure.
figure,plot(s);
hold on,plot(cell2mat(circshift(threshold_cfar,G)),'r--','LineWidth',2)
%hold on, plot (cell2mat(circshift(signal_cfar,(T+G))),'g--','LineWidth',4);
hold on, plot (cell2mat(circshift(signal_cfar,(0))),'g--','LineWidth',4);
legend('Signal','CFAR Threshold','detection')