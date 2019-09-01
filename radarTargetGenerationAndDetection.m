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
v=-20e2; %velocity of the target +-> approaching
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
Nd=2;                   % #of doppler cells OR #of sent periods % number of chirps
%The number of samples on each chirp.
Nr=2^(25-Nd);                  %for length of time OR # of range cells
Nr=2^4;                  %for length of time OR # of range cells
% Timestamp for running the displacement scenario for every sample on each chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples
deltaT=t(2)-t(1);
%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal
%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some tests
% https://www.gaussianwaves.com/2014/07/chirp-signal-frequency-sweeping-fft-and-power-spectral-density
% TxRef = chirp(t,fc,Tchirp,fc+B);
% Tx = cos((fc+slope*t/2).*t*2*pi);
% norm(Tx-TxRef)
n=0;
endFreqTx=fc+B;
for i=1:length(t)
    tScal=t(i)-(n*Tchirp);
    currentFreqTx=fc+slope*tScal; % This is 1st derivative of Tx=(fc+slope*t(i)/2)*t(i)*2*pi
    if(currentFreqTx>endFreqTx)
        n=n+1; %next bin
        tScal=t(i)-(n*Tchirp);
        disp(i)
    end   
    Tx(i) = cos((fc+slope*tScal/2)*tScal*2*pi);
    ttilde=tScal-2*(R+v*t(i))/c;
    if(ttilde>0 || n>0)
        Rx(i) = cos((fc*ttilde+slope*ttilde*ttilde/2)*2*pi);
    end
end
n=0;
for i=1:length(t)
    tScal=t(i)-(n*Tchirp);
    currentFreqTx=fc+slope*tScal; % This is 1st derivative of Tx=(fc+slope*t(i)/2)*t(i)*2*pi
    Tx2(i) = cos((fc+slope*tScal/2)*tScal*2*pi);
    ttilde=tScal-2*(R+v*t(i))/c;
    if(ttilde>0 || n>0)
        Rx(i) = cos((fc*ttilde+slope*ttilde*ttilde/2)*2*pi);
    end
    if (mod(i,(Nr))==0)
        n=n+1; %next bin
        disp(i)
    end
end
norm(Tx-Tx2)

Mix=Tx.*Rx;
windowSize=Nr/2^4;
[~,~,~,pxx1,fc1,tc1] = spectrogram(Tx,windowSize,ceil(windowSize*0.7),windowSize,1/deltaT,'yaxis','MinThreshold',-80);
ylim([77,77.2]);
title('Linear Chirp Tx');
windowSize=Nr/2^4;
[~,~,~,pxx2,fc2,tc2] = spectrogram(Rx,windowSize,ceil(windowSize*0.7),windowSize,1/deltaT,'yaxis','MinThreshold',-80);
title('Linear Chirp Rx');
ylim([77,78]);
plot(tc1(pxx1>0),fc1(pxx1>0),'.');
hold on;
plot(tc2(pxx2>0),fc2(pxx2>0),'.r');
legend('Tx','Rx')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Signal generation and Moving Target simulation
% Running the radar scenario over the time.

for i=1:length(t)
    
    
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity.
    range=R+v*t(i);
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal.
    
    
    Tx(i) = cos((fc*t(i)+slope*t(i)*t(i)/2)*2*pi);
    %   tau = (range/c)*2; %trip time
    %   Rx(i) = cos((fc* (t(i)-tau) +(slope*(t(i)-tau)*(t(i)-tau))/2)*2*pi);
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    %  Mix(i) = cos(2*pi*((2*slope*range/c)*t(i)+(2*fc*v/c)*Tchirp));
    
end
% Mix=Tx.*Rx;
%% END Signal generation and Moving Target simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RANGE MEASUREMENT


% *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
S=reshape(Mix,Nr,Nd);

% *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.

% *%TODO* :
% Take the absolute value of FFT output

% *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.


%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1)

% *%TODO* :
% plot FFT output


axis ([0 200 0 1]);



%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

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
noise_level = zeros(1,1);


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



