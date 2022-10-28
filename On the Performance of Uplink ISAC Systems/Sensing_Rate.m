%% ----------------------------------------------------------------------------- %%
% This Matlab codes characterize the sensing rate in the paper:
% Title: On the Performance of Downlink MIMO-ISAC
% Authors: Chongjun Ouyang, Yuanwei Liu, and Hongwen Yang
% Download article: https://arxiv.org/pdf/2209.01028.pdf
%
% License: If you in any way use this code for research that results in publications, 
% please cite our paper as described above.
%% ----------------------------------------------------------------------------- %%
close all;
clear;

%% Simulation Parameters of the Considered ISAC Systems
M = 4;% Number of Transmit Antennas
N = 5;% Number of Receive Antennas
L = 30;% Length of Radar Waveform
K = 4;% Number of Communication Users
kappa = 0.5;% Power allocation factor in FDSAC
mu = 0.5;% Spectrum allocation factor in FDSAC
eig = [1;0.1;0.05;0.01];% Eigenvalues of the correlation matrix of the target response matrix
R = diag(eig);
Monte = 1e5;
%% Generate communication channels
if K-M==0
    h1 = (abs(1/sqrt(2)*randn(K-M+1,Monte)+1j*1/sqrt(2)*randn(K-M+1,Monte))).^2;
    h2 = (abs(1/sqrt(2)*randn(K-M+1,Monte)+1j*1/sqrt(2)*randn(K-M+1,Monte))).^2;
    h3 = (abs(1/sqrt(2)*randn(K-M+1,Monte)+1j*1/sqrt(2)*randn(K-M+1,Monte))).^2;
    h4 = (abs(1/sqrt(2)*randn(K-M+1,Monte)+1j*1/sqrt(2)*randn(K-M+1,Monte))).^2;
else
    h1 = sum((abs(1/sqrt(2)*randn(K-M+1,Monte)+1j*1/sqrt(2)*randn(K-M+1,Monte))).^2);
    h2 = sum((abs(1/sqrt(2)*randn(K-M+1,Monte)+1j*1/sqrt(2)*randn(K-M+1,Monte))).^2);
    h3 = sum((abs(1/sqrt(2)*randn(K-M+1,Monte)+1j*1/sqrt(2)*randn(K-M+1,Monte))).^2);
    h4 = sum((abs(1/sqrt(2)*randn(K-M+1,Monte)+1j*1/sqrt(2)*randn(K-M+1,Monte))).^2);
end
H = transpose([h1;h2;h3;h4]);
SNR_dB = [-15:1:35];% Sensing SNR
SR_SC = ones(1,length(SNR_dB));
SR_CC = ones(1,length(SNR_dB));
SR_EPA = ones(1,length(SNR_dB));
SR_FDSAC = ones(1,length(SNR_dB));
%% Calculate the sensing rate
for index = [1:1:length(SNR_dB)]
    index
    ps = 10^(SNR_dB(index)/10);
    %% Sensing-Centric
    s = waterfill(ps,1./eig'/L);
    SR_SC(index) = N/L*sum(log2(1+s.*(eig')*L));
    %% Equal Power Allocation
    s = ps/M*ones(1,M);
    SR_EPA(index) = N/L*sum(log2(1+s.*(eig')*L));
    %% Communication-Centric
    s = waterfill(ps*ones(Monte,1),1./H);
    SR_CC(index) = mean(real(N/L*sum(transpose(log2(1+s.*(ones(Monte,1)*(eig'))*L)))));
    %% FDSAC
    s = waterfill(ps*(1-mu),1./eig'/L/(1-kappa));
    SR_FDSAC(index) = N/L*(1-kappa)*sum(log2(1+s.*(eig')*L/(1-kappa)));
end
plot(SNR_dB,SR_SC,'-sk');
hold on;
plot(SNR_dB,SR_CC,'-ok');
hold on;
plot(SNR_dB,SR_EPA,'-*k');
hold on;
plot(SNR_dB,SR_FDSAC,'-vk');
hold on;
%% Calculate the asymptotic sensing rate
SR_ISAC_Asym = N*M/L*(log2(10.^(SNR_dB/10))+1/M*sum(log2(L/M*eig)));
plot(SNR_dB,SR_ISAC_Asym,'--k');
hold on;
SR_FDSAC_Asym = N*M/L*(1-kappa)*(log2(10.^(SNR_dB/10))+1/M*sum(log2((1-mu)*L/(1-kappa)/M*eig)));
plot(SNR_dB,SR_FDSAC_Asym,'--k');
ylim([0,7.5]);
xlim([-15,35]);
grid on;
Data = [SR_SC;SR_CC;SR_EPA;SR_FDSAC;SR_ISAC_Asym;SR_FDSAC_Asym];
save Sensing_Rate_Data Data