%% ----------------------------------------------------------------------------- %%
% This Matlab codes characterize the uplink sensing rate in the paper:
% Title: Performance of Downlink and Uplink Integrated Sensing and Communications (ISAC) Systems
% Authors: Chongjun Ouyang, Yuanwei Liu, and Hongwen Yang
% IEEE Xplore: https://ieeexplore.ieee.org/document/9800940
% Download article: https://arxiv.org/pdf/2202.06207v2.pdf
%
% License: If you in any way use this code for research that results in publications, 
% please cite our paper as described above.
%% ----------------------------------------------------------------------------- %%
close all;
clear;

%% Simulation Parameters of the Considered ISAC Systems
N = 2; % Number of Receiver Antennas
M = 2; % Number of Transmit Antennas
L = 4; % Length of Radar Waveform
K = 2; % Number of Communication Users
Rt = ones(M,M); % Correlation matrix of the target response
for i = [1:1:M]
    for j = [1:1:M]
        Rt(i,j)=0.7^(abs(i-j));
    end
end
eig_t = svd(Rt);
alpha = 0.5; % Spectrum allocation factor in FDSAC
SR = ones(1,length(SNR_dB));
SR1 = ones(1,length(SNR_dB));
for index = [1:1:length(SNR_dB)]
    ps = 10^(SNR_dB(index)/10);
    s = waterfill(ps,1./eig_t');
    SR(index) = N/L*sum(log2(1+s.*(eig_t'))); % Sensing rate in ISAC
    ps = 10^(SNR_dB(index)/10)/alpha;
    s = waterfill(ps,1./eig_t');
    SR1(index) = alpha*N/L*sum(log2(1+s.*(eig_t'))); % Sensing rate in FDSAC
end
plot(SNR_dB,SR,'-o');
hold on;
plot(SNR_dB,SR1,'-s');
hold on;
plot(SNR_dB,N/L*sum(log2(eig_t/M))+N*M/L*log2(10.^(SNR_dB/10)),'--');
hold on;
plot(SNR_dB,alpha*(N/L*sum(log2(eig_t/M))+N*M/L*log2(10.^(SNR_dB/10)/alpha)),'--');
ylim([0,15]);
Data = [SR;SR1;N/L*sum(log2(eig_t/M))+N*M/L*log2(10.^(SNR_dB/10));alpha*(N/L*sum(log2(eig_t/M))+N*M/L*log2(10.^(SNR_dB/10)/alpha))];
save Uplink_Sensing_Rate Data