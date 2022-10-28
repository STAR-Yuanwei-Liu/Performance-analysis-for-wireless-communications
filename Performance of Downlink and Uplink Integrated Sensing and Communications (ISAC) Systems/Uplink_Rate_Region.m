%% ----------------------------------------------------------------------------- %%
% This Matlab codes characterize the uplink rate region in the paper:
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
Pc = 5; % Communication SNR [dB]
Ps = 10; % Sensing SNR [dB]
resolution = [0:0.001:1];
FDSAC_ECR = zeros(1,length(resolution));
FDSAC_SR = zeros(1,length(resolution));
ISAC_ECR = zeros(1,length(resolution));
ISAC_SR = zeros(1,length(resolution));
Monte_Carlo = 1e4;
%% Calculate the sensing rate of FDSAC
%% Calculate the sensing rate of ISAC
for index = [1:1:length(resolution)]
    alpha = resolution(index);
    ps = 10^(Ps/10);
    pc = 10^(Pc/10);
    if alpha~=1
        s = waterfill(ps,(1-alpha)./(eig_t'));
        FDSAC_SR(index) = (1-alpha)*N/L*sum(log2(1+1/(1-alpha)*s.*(eig_t'))); % Sensing rate in FDSAC
    end
    alpha = resolution(index);
    ps = alpha*10^(Ps/10);
    pc = 10^(Pc/10);
    s = waterfill(ps,1./(eig_t'));
    ISAC_SR(index) = N/L*sum(log2(1+s.*(eig_t'))); % Sensing rate in ISAC
end
pc = 10^(Pc/10);
%% Calculate the communication rate in ISAC and FDSAC
parfor C_index = [1:1:Monte_Carlo]
    H = 1/sqrt(2)*randn(K,N) + 1j*1/sqrt(2)*randn(K,N);
    FDSAC_ECR_Tmp = zeros(1,length(resolution));
    ISAC_ECR_Tmp = zeros(1,length(resolution));
    for index = [1:1:length(resolution)]
        [index,C_index]
        alpha = resolution(index);
        if alpha~=0
            FDSAC_ECR_Tmp(index) = alpha*real(log2(det(eye(K)+pc/alpha*(H*H')))); % Communication rate in FDSAC
        end
        ps = alpha*10^(Ps/10);
        s = waterfill(ps,1./(eig_t'));
        Radar_I = zeros(L,1);
        Radar_I([1:1:M]) = s'.*eig_t;
        sigma = 1 + Radar_I;
        % Communication rate in ISAC
        ISAC_ECR_Tmp(index) = (real(log2(det(eye(K)+pc/sigma(1)*(H*H')))) + real(log2(det(eye(K)+pc/sigma(2)*(H*H')))) +...
                              real(log2(det(eye(K)+pc/sigma(3)*(H*H')))) + real(log2(det(eye(K)+pc/sigma(4)*(H*H')))))/L;
    end
    FDSAC_ECR = FDSAC_ECR + FDSAC_ECR_Tmp;
    ISAC_ECR = ISAC_ECR + ISAC_ECR_Tmp;
end
FDSAC_ECR = FDSAC_ECR/Monte_Carlo; 
ISAC_ECR = ISAC_ECR/Monte_Carlo;
plot(ISAC_ECR,ISAC_SR,':');
hold on;
plot(FDSAC_ECR,FDSAC_SR,'-.k');
hold on;
plot([0,ISAC_ECR(end)],[ISAC_SR(end),ISAC_SR(end)],'k-');
hold on;
plot([ISAC_ECR(end),ISAC_ECR(end)],[0,ISAC_SR(end)],'k-');
hold on;
Data = [ISAC_ECR;ISAC_SR;FDSAC_ECR;FDSAC_SR];
save Uplink_Rate_Region Data