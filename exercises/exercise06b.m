%% Solution by Stefan Sicklinger 2019: CA-CFAR
clc;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data_points
Ni = 50;
Nj = 50;
% Generate random noise
s=abs(rand(Ni,Nj));
%Targets location. Assigning bin 100, 200, 300 and 700 as Targets
%with the amplitudes of 8, 9, 4, 11.
s(ceil(0.35*Ni),ceil(0.35*Nj))=8;
s(ceil(0.25*Ni),ceil(0.65*Nj))=9;
s(ceil(0.3*Ni),ceil(0.3*Nj))=5;
s(ceil(0.7*Ni),ceil(0.7*Nj))=11;
%Training Cells
T=[8,4];
%Guard Cells
G=[4,2];
% 2. Slide window across the signal length
% initlal image
CFARImage=zeros(Ni,Nj);
threshold=zeros(Ni,Nj);
singalFiltered=zeros(Ni,Nj);
% Offset : Adding room above noise threshold for desired SNR
offset=10;
% 3 -> training cell
% 2 -> guard cell
% 1 -> cell under test
iSpan=[T(1)+G(1), Ni-T(1)-G(1)+1 ];
jSpan=[T(2)+G(2), Nj-T(2)-G(2)+1 ];
for j = 1:Nj
    for i = 1:Ni
        if(~(i>iSpan(1) && i<iSpan(2) && j>jSpan(1) && j<jSpan(2)))
            s(i,j)=0;
            CFARImage(i,j)=-1;
        end
    end
end
imagesc(CFARImage);
colorbar;
caxis([-1 3]);
for j = 1:Nj
    for i = 1:Ni
        noise_level=0;
        lowerIndexi=i-G(1)-T(1);
        lowerIndexj=j-G(2)-T(2);
        upperIndexi=i+G(1)+T(1);
        upperIndexj=j+G(2)+T(2);
        CFARImage(i,j)=1; % CUT
        if((lowerIndexi>0 && lowerIndexj>0) && (upperIndexi<=Ni && upperIndexj<=Nj) )
            for k = lowerIndexi:upperIndexi
                for l = lowerIndexj:upperIndexj
                    if(((i~=k) || (j~=l)))% spare CUT
                        CFARImage(k,l)=3; % coloring training cell
                        noise_level=s(k,l)+noise_level;
                        addSignal=addSignal+1;
                    end
                end
            end
            % remove guard
            for k = (i-G(1)):(i+G(1))
                for l = (j-G(2)):(j+G(2))
                    if(((i~=k) || (j~=l)))% spare CUT
                        CFARImage(k,l)=CFARImage(k,l)-1; % coloring guard cell
                        noise_level=noise_level-s(k,l);
                        minusSignal=minusSignal+1;
                    end
                end
            end
        end
        %Compute CFAR
        totalNumberTrainingCells=(2*(T(1)+G(1))+1)*(2*(T(2)+G(2))+1)-((2*G(1)+1)*(2*G(2)+1));
        threshold(i,j) = (noise_level/(totalNumberTrainingCells))*offset;
        %Filter the signal above the threshold
        signal =s(i,j);
        if(s(i,j)<threshold(i,j))
            signal=0;
        end
        singalFiltered(i,j) = signal;
        imagesc(CFARImage);
        colorbar;
        caxis([-1 3]);
        pause(0);
        CFARImage=zeros(Ni,Nj);
        for m = 1:Nj
            for n = 1:Ni
                if(s(n,m)==0)
                    CFARImage(n,m)=-1;
                end
            end
        end
    end
end
%plot the output
subplot(3,1,1);
surf(s);
title('Raw signal');
subplot(3,1,2);
surf(threshold);
title('CA CFAR threshold');
subplot(3,1,3);
surf(singalFiltered);
title('CA CFAR filtered');