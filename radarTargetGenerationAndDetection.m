clear all;
format long;
clc;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Radar Specifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
fc= 77e9;
% Max Range = 200m
Rmax= 200;
% Range Resolution = 1 m
dres=1;
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%speed of light = 3e8
%% User Defined Range and Velocity of target
% define the target's initial position and velocity. Note : Velocity
% remains contant
R=110; %initial distance of the target
v=-20; %velocity of the target - -> approaching deltaFD = positiv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FMCW Waveform Generation
%Speed of light (m/s)
c = 3*10^8;
%Calculate the wavelength (m)
lambda = c/fc;
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.
B=c/(2*dres); %(Hz)
Tchirp=5.5*2*Rmax/c; %should be at least 5 to 6 times the round trip time (s)
slope=B/Tchirp; %(-)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END FMCW Waveform Generation
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT for Doppler Estimation.
Nd=2^6;                   % #of doppler cells OR #of sent periods % number of chirps
%The number of samples on each chirp.
Nr=2^(25-Nd); %for length of time OR # of range cells
Nr=2^12;
% Timestamp for running the displacement scenario for every sample on each chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples
deltaT=t(2)-t(1);
%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal
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
endFreqRx=endFreqTx;
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
    ttilde=tScal2-2*(R+v*t(i))/c;
    currentFreqRx=fc+slope*ttilde;  
    if(currentFreqRx>endFreqRx)
        m=m+1; %next bin
        tScal2=t(i)-(m*Tchirp);
        ttilde=tScal2-2*(R+v*t(i))/c;
    end
    if(ttilde>0)
        Rx(i) = cos((fc+slope*ttilde/2)*ttilde*2*pi);
    end
end
% n=0;
% for i=1:length(t)
%       r_t(i)=R+v*t(i); % range of the target in terms of its velocity and initial range
%       ta(i)=2*r_t(i)/c; % delay for received signal
%       if i>n*Nr && i<=(n+1)*Nr % doing this for D of periods (nt length of pulse)
%           Tx(i)=sin(2*pi*(fc*t(i)+.5*slope*t(i)^2-slope*t(i)*n*Tchirp)); %transmitted signal
%           Rx(i)=sin(2*pi*(fc*(t(i)-ta(i))+.5*slope*(t(i)-ta(i))^2-slope*(t(i)-ta(i))*n*Tchirp)); %received signal
%       else
%           n=n+1;
%           Tx(i)=sin(2*pi*(fc*t(i)+.5*slope*t(i)^2-slope*t(i)*n*Tchirp)); %transmitted signal
%           Rx(i)=sin(2*pi*(fc*(t(i)-ta(i))+.5*slope*(t(i)-ta(i))^2-slope*(t(i)-ta(i))*n*Tchirp)); %received signal
%       end
% end
Mix=Tx.*Rx;
% windowSize=Nr/2^4;
% [~,~,~,pxx1,fc1,tc1] = spectrogram(Tx,windowSize,ceil(windowSize*0.7),windowSize,1/deltaT,'yaxis','MinThreshold',-80);
% ylim([77,77.2]);
% title('Linear Chirp Tx');
% windowSize=Nr/2^4;
% [~,~,~,pxx2,fc2,tc2] = spectrogram(Rx,windowSize,ceil(windowSize*0.7),windowSize,1/deltaT,'yaxis','MinThreshold',-80);
% title('Linear Chirp Rx');
% ylim([77,78]);
% x1=tc1(pxx1>0);
% y1=fc1(pxx1>0);
% plot(x1,y1,'.');
% hold on;
% x2=tc2(pxx2>0);
% y2=fc2(pxx2>0);
% plot(x2,y2,'.r');
% legend('Tx','Rx');
% disp(['================================================================']);
% disp(['VERIFICATION']);
% disp(['================================================================']);
% disp(['Start freq TX ref: ' num2str(fc) ' (Hz), ' num2str(fc/1e9) ' (GHz)']);
% disp(['Start freq TX: ' num2str(y1(1)) ' (Hz), ' num2str(y1(1)/1e9) ' (GHz) , error (%): ' num2str(((y1(1)-fc)/fc)*1e2)]);
% disp(['End freq TX ref: ' num2str(endFreqTx) ' (Hz), ' num2str(endFreqTx/1e9) ' (GHz)']);
% disp(['================================================================']);
% startFreqRxRef=(c/(c+v))*fc;
% disp(['Start freq RX ref: ' num2str(startFreqRxRef) ' (Hz), ' num2str(startFreqRxRef/1e9)  ' (GHz)']);
% startFreqRx=y2(1);
% disp(['Start freq RX: ' num2str(y2(1)) ' (Hz), ' num2str(startFreqRx/1e9) ' (GHz), error (%): ' num2str(((startFreqRxRef-startFreqRx)/startFreqRx)*1e2)]);
% disp(['End freq RX ref: ' num2str(endFreqRx) ' (Hz), ' num2str(endFreqRx/1e9) ' (GHz)']);
% disp(['================================================================']);
% disp(['Doppler freq shift ref: ' num2str(deltaDoppler) ' (Hz), ' num2str(deltaDoppler/1e9) ' (GHz)']);
% %disp(['Doppler freq shift: ' num2str(y2(1)-y1(2)) ' (Hz), ' num2str((y2(1)-y1(2))/1e9) ' (GHz), error (%): ' num2str((((y2(1)-y1(2))-deltaDoppler)/deltaDoppler)*1e2)]);
%% END Signal generation and Moving Target simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RANGE MEASUREMENT

%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
S=reshape(Mix,Nr,Nd);
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
Y = fft(S(:,1));
% Take the absolute value of FFT output
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
P2 = abs(Y/Nr);
P2 = fftshift(P2);
P1 = P2(Nr/2+1:end);
%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1)
% plot FFT output
dist = linspace(0,(Nr/2)-1,Nr/2);
plot(dist,P1);
xlabel('Distance (m)');
ylabel('|P1(f)|')
axis ([0 Rmax 0 1]);


%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


m1=reshape(Mix,Nr,Nd); %generating matrix ---> each row showing range info for one period AND each column showing number of periods
[My,Ny]=size(m1');
win=hamming(Ny);
m2=conj(m1).*(win*ones(1,My)); %taking conjugate and applying window for sidelobe reduction (in time domain)
Win=fft(hamming(My),Nd);
M2=(fft(m2,2*Nr)); %First FFT for range information
M3=fftshift(fft(M2',2*Nd)); %Second FFT for doppler information
[My,Ny]=size(M3);
doppler=linspace(-Nd,Nd,My);
range=linspace(-Nr,Nr,Ny);
figure;contour(range,doppler,abs(M3));grid on
xlabel('Range')
ylabel('Doppler')
figure;mesh(range,doppler,abs(M3))
xlabel('Range')
ylabel('Doppler')

% Range Doppler Map Generation.
% This code is buggy therefore; I needed to reimpelemnt 
%%%%%%%%%%%%%%%%
Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);

%%%%%%%%%%%%%

Mix=reshape(Mix,Nr,Nd);
% 2D FFT using the FFT size for both dimensions.
P2=(fft(Mix,Nr)); %First FFT for range information
P2=fftshift(fft(P2',Nd)); %Second FFT for doppler information
P1 = P2(:,Nr/2+1:end);
sig_fft2 = fft2(Mix,Nr,Nd);
RDM = abs(P1);
%RDM = 10*log10(P2) ;
dist = linspace(0,(Nr/2)-1,Nr/2);
doppler=linspace(-Nd,Nd,Nd);
figure();
surf(dist,doppler,RDM);
%use the surf function to plot the output of 2DFFT and to show axis in both

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions.

% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under
%test (CUT) for accurate estimation

% *%TODO* :
% offset the threshold by SNR value in dB

% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
mwnoise_level = zeros(1,1);


% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.


% Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
% CFAR





% *%TODO* :
% The process above will generate a thresholded block, which is smaller
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0.









% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis,'replace this with output');
colorbar;



