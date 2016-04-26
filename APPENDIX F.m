% This code is to plot receiver operating characteristic curve  
% Cooperative Spectrum Sensing using Joint Estimation and Detection,
% when the primary signal is real Gaussian signal and noise is
% addive white real Gaussian. Here, the threshold is available
% analytically. Softened Hard(3-Bit Implementation). SNR vs Pd.
% Code written by: Akshay Prasad, Shiv Nadar University India.
% Refer to Three Bits Softened Decision Scheme in Cooperative Spectrum Sensing Among
% Cognitive Radio Networks by Hefdhallah Sakran et al.
clc
close all
clear all
L = 1000;
no = 10; %no. of users to take
snr_dB =  -16:0.5:-4;  % SNR in decibels
alpha =[.84 .81 .76 .68 .54 .13];
g = [1 .33 .42 .6 .75 1.5 3];
Pf=zeros(1,7);
Pf(1)=0.01;
for cnti = 1:6
    Pf(cnti+1)=Pf(cnti)*alpha(cnti);
end

%Creating 15 signal streams being recieved     
y=zeros(15,1000);
n=zeros(15,1000);
s=zeros(15,1000);
%Adding Noise to the 15 signal streams
dec=zeros(1,no);%Array to store the decisions of the 10 independednt sensors
energy=zeros(no,1000);%Calculating Energy
energy_fin=zeros(no,1);
thresh=zeros(7,1);
for m = 1:length(snr_dB)
    i = 0;
    snr = 10.^(snr_dB(m)./10); % Linear Value of SNR
for ll = 1:15 
   n(ll,:) = randn(1,L); %AWGN noise with mean 0 and variance 1
   s(ll,:) = sqrt(snr).*randn(1,L); % Real valued Gaussian Primary User Signal 
   y(ll,:) = s(ll,:) + n(ll,:) +1; % Received signal at SU
end;
correl=zeros(15,15); %Calculating the correlation between the sensors
for mm = 1:15 %To Calculate Mutual Information to Extract Entropy
    for ff = 1:15
        correl(mm,ff)=mutualInformation(y(mm,:),y(ff,:));
    end
end;
endarr=findminimum(correl,15);%Finds the index of the sensors that have minimum correlation
 for kk=1:1000 % Number of Monte Carlo Simulations
  for ll = 1:15 
   n(ll,:) = randn(1,L); %AWGN noise with mean 0 and variance 1
   s(ll,:) = sqrt(snr).*randn(1,L); % Real valued Gaussian Primary User Signal 
   y(ll,:) = s(ll,:) + n(ll,:); % Received signal at SU
  end
  for ff=1:no
   energy(ff,:) = abs(y((endarr(ff)),:)).^2; % Energy of received signals over N samples
   energy_fin(ff,1) =(1/L).*sum(energy(ff,:)); % Test Statistic for the energy detection
  end
  for xx=1:7
      thresh(xx)= (qfuncinv(Pf(xx))./sqrt(L))+1; 
  end
  ni=zeros(1,7);
  for count =1:no
      if(thresh(1)<energy_fin(count,1)&&thresh(2)>energy_fin(count,1))
          ni(1)=ni(1)+1;
      elseif(thresh(2)<energy_fin(count,1)&&thresh(3)>energy_fin(count,1))
          ni(2)=ni(2)+1;
      elseif(thresh(3)<energy_fin(count,1)&&thresh(4)>energy_fin(count,1))
          ni(3)=ni(3)+1;
      elseif(thresh(4)<energy_fin(count,1)&&thresh(5)>energy_fin(count,1))
          ni(4)=ni(4)+1;
      elseif(thresh(5)<energy_fin(count,1)&&thresh(6)>energy_fin(count,1))
          ni(5)=ni(5)+1;
      elseif(thresh(6)<energy_fin(count,1)&&thresh(7)>energy_fin(count,1))
          ni(6)=ni(6)+1;
      elseif(thresh(7)<energy_fin(count,1))
          ni(7)=ni(7)+1; 
      end      
  end
  test=ni(1)*g(1)+ni(2)*g(2)+ni(3)*g(3)+ni(4)*g(4)+ni(5)*g(5)+ni(6)*g(6)+ni(7)*g(7);
  
  if(test>3)
      i=i+1;
  end
  
 end
disp(i);
Pd(m) = i/kk; 
end
plot(snr_dB, Pd)
hold on
title('SNR vs Pd(Softened - Hard Fusion Rule)');
xlabel('SNR') % x-axis label
ylabel('Pd (Probablity of Detection)') % y-axis label
