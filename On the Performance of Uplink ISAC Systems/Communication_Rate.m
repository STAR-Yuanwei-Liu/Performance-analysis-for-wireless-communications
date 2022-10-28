%% ----------------------------------------------------------------------------- %%
% This Matlab codes characterize the ergodic communication rate in the paper:
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
SNR_dB = [-5:1:35];% Communication SNR
CR_SC = ones(1,length(SNR_dB));
CR_CC = ones(1,length(SNR_dB));
CR_EPA = ones(1,length(SNR_dB));
CR_FDSAC = ones(1,length(SNR_dB));
%% Calculate the communication rate
for index = [1:1:length(SNR_dB)]
    ps = 10^(SNR_dB(index)/10);
    %% Sensing-Centric
    s = waterfill(ps,1./eig'/L);
    CR_SC(index) = mean(sum(transpose(log2(1+(ones(Monte,1)*s).*H))));
    %% Equal Power Allocation
    s = ps/M*ones(1,M);
    CR_EPA(index) = mean(sum(transpose(log2(1+(ones(Monte,1)*s).*H))));
    %% Communication-Centric
    s = waterfill(ps*ones(Monte,1),1./H);
    CR_CC(index) = mean(sum(transpose(log2(1+s.*H))));
    %% FDSAC
    s = waterfill(ps*mu,1./H/kappa);
    CR_FDSAC(index) = kappa*mean(sum(transpose(log2(1+s.*H/kappa))));
end
plot(SNR_dB,CR_SC,'-sk');
hold on;
plot(SNR_dB,CR_CC,'-ok');
hold on;
plot(SNR_dB,CR_EPA,'-*k');
hold on;
plot(SNR_dB,CR_FDSAC,'-vk');
hold on;
%% Calculate the asymptotic ergodic communication rate
CR_ISAC_Asym = M*(log2(10.^(SNR_dB/10))-log2(M)+psi(K-M+1)/log(2));
plot(SNR_dB,CR_ISAC_Asym,'--k');
hold on;
CR_FDSAC_Asym = (1-kappa)*M*(log2(10.^(SNR_dB/10))-log2(M)+psi(K-M+1)/log(2));
plot(SNR_dB,CR_FDSAC_Asym,'--k');
ylim([0,35.5]);
xlim([-5,35]);
grid on;
Data = [CR_SC;CR_CC;CR_EPA;CR_FDSAC;CR_ISAC_Asym;CR_FDSAC_Asym];
save Communication_Rate_Data Data