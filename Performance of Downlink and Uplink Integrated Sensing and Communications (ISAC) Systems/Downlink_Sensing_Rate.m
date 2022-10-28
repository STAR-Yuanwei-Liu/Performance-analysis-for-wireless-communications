%% ----------------------------------------------------------------------------- %%
% This Matlab codes characterize the downlink sensing rate in the paper:
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
Rt = ones(M,M); % Correlation matrix for target response
for i = [1:1:M]
    for j = [1:1:M]
        Rt(i,j)=0.7^(abs(i-j));
    end
end
R = ones(M,M); % Correlation matrix for communication user
for i = [1:1:M]
    for j = [1:1:M]
        R(i,j)=0.8^(abs(i-j));
    end
end
eig_t = svd(Rt);
Pc = 5;
pc = 10^(Pc/10);
p_dpc = pc/K*[1,1];
Monte_Carlo = 5e4;
DPC_Cov = 0;
%% Obtain the convariance matrix for communications
for index = [1:1:Monte_Carlo]
    index
    H = (1/sqrt(2)*randn(K,M) + 1j*1/sqrt(2)*randn(K,M))*(R^0.5); % Correlated Rayleigh Fading
    h1 = H(1,:); % Channel of communication user 1
    h2 = H(2,:); % Channel of communication user 2
    v = real([h1*h1',h2*h2']);
    iteration = 40;
    %% Iterative Water-Filling Method
    % For more details, please refer to the following paper
    % Nihar Jindal, Wonjong Rhee, Sriram Vishwanath, Syed Ali Jafar, and Andrea Goldsmith, 
    % "Sum Power Iterative Water-filling for Multi-Antenna Gaussian Broadcast Channels," 
    % IEEE Transactions on Information Theory, Vol. 51, No. 4, pp. 1570-1580, April 2005.
    for j = [1:1:iteration]    
        p_dpc1 = waterfill(pc,1./(v - (p_dpc*((abs(h1*h2'))^2))./(1+p_dpc.*wrev(v))));
        p_dpc = 1/2*(p_dpc1+p_dpc);
    end
    B1 = eye(M) - h2'*h2*p_dpc(2)/(1+p_dpc(2)*h2*h2');
    M1 = B1*(h1'*h1)*B1*p_dpc(1)/(h1*B1*h1');
    M2 = (h2'*h2)*p_dpc(2)*(1+h2*M1*h2')/(h2*h2');
    DPC_Cov = DPC_Cov + M1 + M2;
end
DPC_Cov = DPC_Cov/Monte_Carlo;
%% Calculate the sensing rate
SNR_dB = [-10:5:40];
alpha = 0.5;
SR = ones(1,length(SNR_dB));
SR1 = ones(1,length(SNR_dB));
for index = [1:1:length(SNR_dB)]
    ps = 10^(SNR_dB(index)/10);
    s = waterfill(ps,(1 + real(trace(Rt*DPC_Cov)))./eig_t');
    SR(index) = N/L*sum(log2(1+1/(1 + real(trace(Rt*DPC_Cov)))*s.*(eig_t'))); %Sensing rate for ISAC
    ps = 10^(SNR_dB(index)/10)/alpha;
    s = waterfill(ps,1./eig_t');
    SR1(index) = alpha*N/L*sum(log2(1+s.*(eig_t'))); %Sensing rate for FDSAC
end
plot(SNR_dB,SR,'-o');
hold on;
plot(SNR_dB,SR1,'-s');
hold on;
plot(SNR_dB,N/L*sum(log2(eig_t/M/(1 + real(trace(Rt*DPC_Cov)))))+N*M/L*log2(10.^(SNR_dB/10)),'--');
hold on;
plot(SNR_dB,alpha*(N/L*sum(log2(eig_t/M))+N*M/L*log2(10.^(SNR_dB/10)/alpha)),'--');
ylim([0,15]);
Data = [SR;SR1;N/L*sum(log2(eig_t/M/(1 + real(trace(Rt*DPC_Cov)))))+N*M/L*log2(10.^(SNR_dB/10));alpha*(N/L*sum(log2(eig_t/M))+N*M/L*log2(10.^(SNR_dB/10)/alpha))];
save Downlink_Sensing_Rate Data