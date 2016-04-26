% This code is to plot SNR vs Pd characteristic curve  
% Cooperative Spectrum Sensing using Joint Estimation and Detection,
% when the primary signal is real Gaussian signal and noise is
% addive white real Gaussian. Here, the threshold is available
% analytically and soft fusion is carried out using EGC.
% Code written by: Akshay Prasad, Shiv Nadar University India.
clc
close all
clear all
M=1;
L = 1000;
no = 4; %no. of users to take
snr_dB =  -16:0.5:-4;  % SNR in decibels
Pf = 0.01; % Pf = Probability of False Alarm
%Creating 15 signal streams being recieved     
y=zeros(15,1000);
n=zeros(15,1000);
s=zeros(15,1000);
%Adding Noise to the 15 signal streams
dec=zeros(1,no);%Array to store the decisions of the 10 independednt sensors
energy=zeros(no,1000);%Calculating Energy
energy_fin=0;
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
  avg=zeros(1,1000);
  for ff=1:no
   energy(ff,:) = abs(y((endarr(ff)),:)).^2; % Energy of received signals over N samples
   avg(1,:)=avg(1,:)+energy(ff,:);
  end
  avg=avg./sqrt(no);
  energy_fin =(1/L).*sum(avg(1,:)); % Test Statistic for the energy detection
  thresh = ((qfuncinv(Pf))*sqrt(2*M))+M*sqrt(1/no); % Theoretical value of Threshold, refer, (EGC)Soft Combination and Detection for Cooperative Spectrum Sensing in Cognitive Radio Networks Jun Ma and Ye (Geoffrey) Li
  if(energy_fin>=thresh)
      i=i+1;
  end
 end
disp(i);
Pd(m) = i/kk; 
end
plot(snr_dB, Pd)
hold on
title('SNR vs Pd(Soft Combination Rule)');
xlabel('SNR') % x-axis label
ylabel('Pd (Probablity of Detection)') % y-axis label
