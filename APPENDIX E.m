% This code is to plot receiver operating characteristic curve  
% Cooperative Spectrum Sensing using Joint Estimation and Detection,
% when the primary signal is real Gaussian signal and no(abc)ise is
% addive white real Gaussian. Here, the threshold is available
% analytically.  SNR vs Pd for varying No. of correlated users.
% Code written by: Akshay Prasad, Shiv Nadar University India.
clc
close all
clear all
L = 1000;
snr_dB =  -16:0.5:-4;  % SNR in decibels
Pf = 0.01; % Pf = Probability of False Alarm
%Creating 50 signal streams being recieved     
y=zeros(50,1000);
n=zeros(50,1000);
s=zeros(50,1000);
no=10:1:32;
abc=1;
for abc = 1:length(no)
%Adding no(abc)ise to the 50 signal streams
dec=zeros(1,no(abc));%Array to store the decisions of the 10 independednt sensors
energy=zeros(no(abc),1000);%Calculating Energy
energy_fin=zeros(no(abc),1); 
mag=zeros(1,length(snr_dB));
for m = 1:length(snr_dB)
    i = 0;
    snr = 10.^(snr_dB(m)./10); % Linear Value of SNR
for ll = 1:50 
   n(ll,:) = randn(1,L); %AWGN no(abc)ise with mean 0 and variance 1
   s(ll,:) = sqrt(snr).*randn(1,L); % Real valued Gaussian Primary User Signal 
   y(ll,:) = s(ll,:) + n(ll,:) +1; % Received signal at SU
end;
correl=zeros(50,50); %Calculating the correlation between the sensors
for mm = 1:50 %To Calculate Mutual Information to Extract Entropy
    for ff = 1:50
        correl(mm,ff)=mutualInformation(y(mm,:),y(ff,:));
    end
end;
for ctrn=1:2:no(abc)-1
    mag(1,m)=mag(1,m)+mutualInformation(y(ctrn,:),y(ctrn+1,:));
end
endarr=findminimum(correl,50);%Finds the index of the sensors that have minimum correlation
 for kk=1:1000 % Number of Monte Carlo Simulations
  for ll = 1:50 
   n(ll,:) = randn(1,L); %AWGN no(abc)ise with mean 0 and variance 1
   s(ll,:) = sqrt(snr).*randn(1,L); % Real valued Gaussian Primary User Signal 
   y(ll,:) = s(ll,:) + n(ll,:); % Received signal at SU
  end
  for ff=1:no(abc)
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
disp(mean(mag));
subplot(6,4,abc) % first subplot
plot(snr_dB, Pd);
hold on
title(['SNR vs Pd, L=' num2str(no(abc)) 'Correl='  num2str(mean(mag))]);
xlabel('SNR') % x-axis label
ylabel('Pd') % y-axis label
end