%% ----------------------------------------------------------------------------- %%
% This Matlab codes characterize the outage probability of communication rate in the paper:
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
SNR_dB = [-5:5:35];% Communication SNR
Monte1 = 1e4;
OP_SC = zeros(1,length(SNR_dB));
OP_CC = zeros(1,length(SNR_dB));
OP_EPA = zeros(1,length(SNR_dB));
OP_FDSAC = zeros(1,length(SNR_dB));
Monte = 1e4;
for outer_index = [1:1:Monte1]
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
    OP_SC_tmp = zeros(1,length(SNR_dB));
    OP_CC_tmp = zeros(1,length(SNR_dB));
    OP_EPA_tmp = zeros(1,length(SNR_dB));
    OP_FDSAC_tmp = zeros(1,length(SNR_dB));
    %% Calculate the outage probability
    for index = [1:1:length(SNR_dB)]
        [outer_index,index]
        ps = 10^(SNR_dB(index)/10);
        %% Sensing-Centric
        s = waterfill(ps,1./eig'/L);
        OP_SC_tmp(index) = length(find((sum(transpose(log2(1+(ones(Monte,1)*s).*H))))<2));
        %% Equal Power Allocation
        s = ps/M*ones(1,M);
        OP_EPA_tmp(index) = length(find((sum(transpose(log2(1+(ones(Monte,1)*s).*H))))<2));
        %% Communication-Centric
        s = waterfill(ps*ones(Monte,1),1./H);
        OP_CC_tmp(index) = length(find((sum(transpose(log2(1+s.*H))))<2));
        %% FDSAC
        s = waterfill(ps*mu,1./H/kappa);
        OP_FDSAC_tmp(index) = length(find(kappa*(sum(transpose(log2(1+s.*H/kappa))))<2));
    end
    OP_SC = OP_SC + OP_SC_tmp;
    OP_EPA = OP_EPA + OP_EPA_tmp;
    OP_CC = OP_CC + OP_CC_tmp;
    OP_FDSAC = OP_FDSAC + OP_FDSAC_tmp;
end
semilogy(SNR_dB,OP_SC/Monte/Monte1,'-sk');
hold on;
semilogy(SNR_dB,OP_CC/Monte/Monte1,'-ok');
hold on;
semilogy(SNR_dB,OP_EPA/Monte/Monte1,'-*k');
hold on;
semilogy(SNR_dB,OP_FDSAC/Monte/Monte1,'-vk');
hold on;
%% Calculate the asymptotic outage probability
semilogy(SNR_dB,1./(((10.^(SNR_dB/10))).^(M*(K-M+1))),'--');
ylim([10^(-6),10^(-0.01)]);
xlim([-5,35]);
grid on;
Data = [OP_SC/Monte/Monte1;OP_CC/Monte/Monte1;OP_EPA/Monte/Monte1;OP_FDSAC/Monte/Monte1;1./(((10.^(SNR_dB/10))).^(M*(K-M+1)))];
save Outage_Probability_Data Data