% This code is to plot receiver operating characteristic curve  
% Cooperative Spectrum Sensing using Joint Estimation and Detection,
% when the primary signal is real Gaussian signal and noise is
% addive white real Gaussian. Here, the threshold is available
% analytically.
% Code written by: Akshay Prasad, Shiv Nadar University India.
clc
close all
clear all
L = 10000;
no = 10; %no. of users to take
snr_dB =  -16:0.5:-4;  % SNR in decibels
Pf = 0.01; % Pf = Probability of False Alarm
%Creating 15 signal streams being recieved     
y=zeros(15,10000);
n=zeros(15,10000);
s=zeros(15,10000);
%Adding Noise to the 15 signal streams
dec=zeros(1,no);%Array to store the decisions of the 10 independednt sensors
energy=zeros(no,10000);%Calculating Energy
energy_fin=zeros(no,1); 
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
 for kk=1:10000 % Number of Monte Carlo Simulations
  for ll = 1:15 
   n(ll,:) = randn(1,L); %AWGN noise with mean 0 and variance 1
   s(ll,:) = sqrt(snr).*randn(1,L); % Real valued Gaussian Primary User Signal 
   y(ll,:) = s(ll,:) + n(ll,:); % Received signal at SU
  end
  for ff=1:no
   energy(ff,:) = abs(y((endarr(ff)),:)).^2; % Energy of received signals over N samples
   energy_fin(ff,1) =(1/L).*sum(energy(ff,:)); % Test Statistic for the energy detection
   thresh = (qfuncinv(Pf)./sqrt(L))+1; % Theoretical value of Threshold, refer, Sensing Throughput Tradeoff in Cognitive Radio, Y. C. Liang
   if(energy_fin(ff,1) >= thresh)  %Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) counter by 1
       dec(ff)=1;
   else
       dec(ff)=0;
   end
  end
  %disp(dec);
  zz=length(dec(dec>0));
  if(zz>=8)  %Modify this to change the Fusion Rule to OR/AND/Majority
      i=i+1;
  end
 end
disp(i);
Pd(m) = i/kk; 
end
plot(snr_dB, Pd)
hold on
title('SNR vs Pd(MAJORITY Rule Fusion)');
xlabel('SNR') % x-axis label
ylabel('Pd (Probablity of Detection)') % y-axis label
