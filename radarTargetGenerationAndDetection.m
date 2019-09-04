clear all;
format long;
clc;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User defined range and velocity of target user input
% define the target's initial position and velocity. Note : Velocity remains contant
R=130; %initial distance of the target
v=-10; %velocity of the target - -> approaching deltaFD = positiv
showWaveGen = 0; % useful for debugging waveform generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Radar specifications user input
% Frequency of operation = 77GHz
fc= 77e9;
% Max Range = 200m
Rmax= 200;
% Range Resolution = 1 m
dres=1;
% Max Velocity = 100 m/s
vmax=70;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FMCW user input waveform generation
%Speed of light (m/s)
c = 3*10^8;
%Calculate the wavelength (m)
lambda = c/fc;
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW chirp using the requirements above.
B=c/(2*dres); %(Hz)
Tchirp=5.5*2*Rmax/c; %should be at least 5 to 6 times the round trip time (s)
slope=B/Tchirp; %(-)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT for Doppler Estimation.
Nd=2^8;                   % #of doppler cells OR #of sent periods % number of chirps
if(showWaveGen)
    Nd=2^2;
end
%The number of samples on each chirp.
Nr=2^12;
if(showWaveGen)
    Nr=2^21;
end
% Timestamp for running the displacement scenario for every sample on each chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples
deltaT=t(2)-t(1);
%Creating the vectors for Tx, Rx
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Signal generation and Moving Target simulation
% Running the radar scenario over the time.
% https://www.gaussianwaves.com/2014/07/chirp-signal-frequency-sweeping-fft-and-power-spectral-density
% some tests
% TxRef = chirp(t,fc,Tchirp,fc+B);
% Tx = cos((fc+slope*t/2).*t*2*pi);
% norm(Tx-TxRef)
n=0;
m=0;
endFreqTx=fc+B;
deltaDoppler=(c/(c+v))*fc-fc;
endFreqRx=endFreqTx+deltaDoppler;
for i=1:length(t)
    tScal1=t(i)-(n*Tchirp);
    currentFreqTx=fc+slope*tScal1; % This is 1st derivative of Tx=(fc+slope*t(i)/2)*t(i)*2*pi
    if(currentFreqTx>endFreqTx)
        n=n+1; %next bin
        tScal1=t(i)-(n*Tchirp);
    end
    Tx(i) = cos((fc+slope*tScal1/2)*tScal1*2*pi);
    %Rx
    tScal2=t(i)-(m*Tchirp);
    ttilde=tScal2-2*(R-v*t(i))/c;
    currentFreqRx=fc+slope*ttilde;
    if(currentFreqRx>endFreqRx)
        m=m+1; %next bin
        tScal2=t(i)-(m*Tchirp);
        ttilde=tScal2-2*(R-v*t(i))/c;
    end
    if(ttilde>0)
        Rx(i) = cos((fc+slope*ttilde/2)*ttilde*2*pi);
    end
end
Mix=Tx.*Rx;
if(showWaveGen)
    windowSize=Nr/2^4;
    [~,~,~,pxx1,fc1,tc1] = spectrogram(Tx,windowSize,ceil(windowSize*0.7),windowSize,1/deltaT,'yaxis','MinThreshold',-80);
    ylim([77,77.2]);
    title('Linear Chirp Tx');
    windowSize=Nr/2^4;
    [~,~,~,pxx2,fc2,tc2] = spectrogram(Rx,windowSize,ceil(windowSize*0.7),windowSize,1/deltaT,'yaxis','MinThreshold',-80);
    title('Linear Chirp Rx');
    ylim([77,78]);
    x1=tc1(pxx1>0);
    y1=fc1(pxx1>0);
    plot(x1,y1,'.');
    hold on;
    x2=tc2(pxx2>0);
    y2=fc2(pxx2>0);
    plot(x2,y2,'.r');
    legend('Tx','Rx');
end
disp(['================================================================']);
disp(['VERIFICATION']);
disp(['================================================================']);
disp(['Start freq TX ref: ' num2str(fc) ' (Hz), ' num2str(fc/1e9) ' (GHz)']);
if(showWaveGen)
    disp(['Start freq TX: ' num2str(y1(1)) ' (Hz), ' num2str(y1(1)/1e9) ' (GHz) , error (%): ' num2str(((y1(1)-fc)/fc)*1e2)]);
end
disp(['End freq TX ref: ' num2str(endFreqTx) ' (Hz), ' num2str(endFreqTx/1e9) ' (GHz)']);
disp(['================================================================']);
startFreqRxRef=(c/(c+v))*fc;
disp(['Start freq RX ref: ' num2str(startFreqRxRef) ' (Hz), ' num2str(startFreqRxRef/1e9)  ' (GHz)']);
if(showWaveGen)
    startFreqRx=y2(1);
    disp(['Start freq RX: ' num2str(y2(1)) ' (Hz), ' num2str(startFreqRx/1e9) ' (GHz), error (%): ' num2str(((startFreqRxRef-startFreqRx)/startFreqRx)*1e2)]);
end
disp(['End freq RX ref: ' num2str(endFreqRx) ' (Hz), ' num2str(endFreqRx/1e9) ' (GHz)']);
disp(['================================================================']);
disp(['Doppler freq shift ref: ' num2str(deltaDoppler) ' (Hz), ' num2str(deltaDoppler/1e9) ' (GHz)']);
if(showWaveGen)
    disp(['Doppler freq shift: ' num2str(y2(1)-y1(2)) ' (Hz), ' num2str((y2(1)-y1(2))/1e9) ' (GHz), error (%): ' num2str((((y2(1)-y1(2))-deltaDoppler)/deltaDoppler)*1e2)]);
    pause(100);
end
%% END Signal generation and Moving Target simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of Range and Doppler FFT respectively.
Mix=reshape(Mix,Nr,Nd);
%% START 1D FFT RANGE MEASUREMENT
% %run the FFT on the beat signal along the range bins dimension (Nr) and
% %normalize.
% Y = fft(Mix(:,1));
% % Take the absolute value of FFT output
% % Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% % Hence we throw out half of the samples.
% P2 = abs(Y/Nr);
% P2 = fftshift(P2);
% P1 = P2(Nr/2+1:end);
% %plotting the range
% figure ('Name','Range from First FFT')
% % plot FFT output
% dist = linspace(0,(Nr/2)-1,Nr/2);
% plot(dist,P1);
% xlabel('Distance (m)');
% ylabel('|Amplitude|');
% xlim([0 Rmax]);
%% END 1D FFT -> this is inefficent and does not apply propper windowing
figure('Name','FFT1D and FFT2D');
Mix_Hann_Wnd=Mix.*(hann(Nr)*ones(1,Nd));%tensor product with sampled hamm window
P2 = fft(Mix_Hann_Wnd);%normalize amplitude spectrum
P2=fft(P2'); % do second FFT on transpose of two-sided spectrum of 1st FFT
P1 = fftshift(P2');
P1 = P1(Nr/2+1:end,:); % Get one-sided spectrum for distance
dist = linspace(0,(Nr/2)-1,Nr/2);
doppler=linspace(-Nd,Nd,Nd)/(Nd/128);
subplot(3,2,1);
plot(dist,abs(P1)); %superpose all distance bins
xlim([0 Rmax]);
title('Distance with Hamm Wnd from 1st FFT');
xlabel('Distance (m)');
ylabel('|Amplitude|');
P2 = fftshift(P2);
subplot(3,2,3);
plot(doppler,abs(P2));
xlim([-vmax vmax]);
title('Velocity with Hamm Wnd from 2nd FFT');
xlabel('Doppler velocity (m/s)');
ylabel('|Amplitude|');
ax5 = subplot(3,2,5);
contour(doppler,dist,abs(P1));
colormap(ax5,hot(8))
xlim([-vmax vmax]);
ylim([0 Rmax]);
title('Contour plot');
xlabel('Doppler velocity (m/s)');
ylabel('Distance (m)');
grid on;
grid minor;
subplot(3,2,[2,4,6]);
surf(doppler,dist,10*log10(abs(P1)));
colormap default
shading interp
title('Surface plot');
xlabel('Doppler velocity (m/s)');
ylabel('Distance (m)');
zlabel('dBW|Amplitude|');
xlim([-vmax vmax]);
ylim([0 Rmax]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CFAR implementation
% Start User input
%Slide Window through the complete Range Doppler Map select the number of Training Cells in both the dimensions.
T=[8,4];
%Select the number of Guard Cells in both dimensions around the Cell under test (CUT) for accurate estimation
G=[4,2];
% offset the threshold by SNR value in dB
offsetThresholdIndB=3;
% End user input
%Create a vector to store noise_level for each iteration on training cells
Ni = Nr/2;
Nj = Nd;
threshold=zeros(Ni,Nj);
singalFiltered=zeros(Ni,Nj);
%The process above will generate a thresholded block, which is smaller than the Range Doppler Map as the CUT cannot be 
% located at the edges of matrix. Hence,few cells will not be thresholded. To keep the map size same set those values to 0.
% --> Setting singal margin to zero
iSpan=[T(1)+G(1), Ni-T(1)-G(1)+1 ];
jSpan=[T(2)+G(2), Nj-T(2)-G(2)+1 ];
for j = 1:Nj
    for i = 1:Ni
        if(~(i>iSpan(1) && i<iSpan(2) && j>jSpan(1) && j<jSpan(2)))
            P1(i,j)=0;
        end
    end
end
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.
RDM=abs(P1);
for j = 1:Nj
    for i = 1:Ni
        noiseLevel=0;
        lowerIndexi=i-G(1)-T(1);
        lowerIndexj=j-G(2)-T(2);
        upperIndexi=i+G(1)+T(1);
        upperIndexj=j+G(2)+T(2);
        if((lowerIndexi>0 && lowerIndexj>0) && (upperIndexi<=Ni && upperIndexj<=Nj) )
            for k = lowerIndexi:upperIndexi
                for l = lowerIndexj:upperIndexj
                    if(((i~=k) || (j~=l)))% spare CUT
                        % coloring training cell
                        noiseLevel=RDM(k,l)+noiseLevel;
                    end
                end
            end
            % remove guard
            for k = (i-G(1)):(i+G(1))
                for l = (j-G(2)):(j+G(2))
                    if(((i~=k) || (j~=l)))% spare CUT
                        % coloring guard cell
                        noiseLevel=noiseLevel-RDM(k,l);
                    end
                end
            end
        end
        %Compute CFAR
        totalNumberTrainingCells=(2*(T(1)+G(1))+1)*(2*(T(2)+G(2))+1)-((2*G(1)+1)*(2*G(2)+1));
        relNoiseLevel=noiseLevel/totalNumberTrainingCells;
        if(relNoiseLevel>0)
            relNoiseLevelOffset=db2pow(10*log10(relNoiseLevel)+offsetThresholdIndB);
            %disp(relNoiseLevelOffset/relNoiseLevel)
            threshold(i,j) = relNoiseLevelOffset;
        else
            threshold(i,j) = relNoiseLevel;
        end
        %Filter the signal above the threshold
        signal =RDM(i,j);
        if(RDM(i,j)<threshold(i,j))
            signal=0;
        end
        singalFiltered(i,j) = signal;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot the output of CA CFAR
figure('Name','CA CFAR #1');
subplot(2,1,1);
surf(doppler,dist,10*log(RDM));
colormap default
shading interp
title('Unfiltered range doppler signal');
xlabel('Doppler velocity (m/s)');
ylabel('Distance (m)');
zlabel('dBW|Amplitude|');
xlim([-vmax vmax]);
ylim([0 Rmax]);
subplot(2,1,2);
surf(doppler,dist,RDM);
colormap default
shading interp
title('Unfiltered range doppler signal');
xlabel('Doppler velocity (m/s)');
ylabel('Distance (m)');
zlabel('|Amplitude|');
xlim([-vmax vmax]);
ylim([0 Rmax]);
figure('Name','CA CFAR #2');
subplot(2,1,1);
surf(doppler,dist,threshold);
colormap default
shading interp
title('CA CFAR threshold');
xlabel('Doppler velocity (m/s)');
ylabel('Distance (m)');
zlabel('|Amplitude|');
subplot(2,1,2);
surf(doppler,dist,singalFiltered);
colormap default
shading interp
title('CA CFAR filtered');
xlabel('Doppler velocity (m/s)');
ylabel('Distance (m)');
zlabel('|Amplitude|');
[M,I] = max(singalFiltered(:));
[i, j] = ind2sub(size(singalFiltered),I);
hold on;
plot3(doppler(j),dist(i),singalFiltered(i,j),'xr','markersize',10);
text(doppler(j),dist(i),singalFiltered(i,j),['  dist: ' num2str(dist(i)) ' (m) doppler: ' num2str(doppler(j)) ' (m/s)'] );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
