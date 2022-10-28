%% ----------------------------------------------------------------------------- %%
% This Matlab codes characterize the downlink ergodic rate of communication users in the paper:
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
R = ones(M,M); % Correlation Matrix of Communication Users
for i = [1:1:M]
    for j = [1:1:M]
        R(i,j)=0.8^(abs(i-j));
    end
end
Monte_Carlo = 1e4; % This should be set to a very large value to make the curve smooth
CSNR_dB = [-10:5:50]; % SNR
Result_ISAC = ones(Monte_Carlo,length(CSNR_dB));
Result_FDSAC = ones(Monte_Carlo,length(CSNR_dB));
Result_ISAC1 = ones(Monte_Carlo,length(CSNR_dB));
Result_FDSAC1 = ones(Monte_Carlo,length(CSNR_dB));
result = ones(1,Monte_Carlo);
alpha = 0.5; % Spectrum allocation factor in FDSAC
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
    for index_dB = [1:1:length(CSNR_dB)]
        pc = 10^(CSNR_dB(index_dB)/10);
        p_dpc = pc/K*[1,1];
        p_fdsac = pc/K*[1,1];
        for j = [1:1:iteration]    
            p_dpc1 = waterfill(pc,1./(v - (p_dpc*((abs(h1*h2'))^2))./(1+p_dpc.*wrev(v))));
            p_fdsac1 = waterfill(pc,1./(v - (p_fdsac*((abs(h1*h2'))^2))./(1+p_fdsac.*wrev(v)))*alpha);
            p_dpc = 1/2*(p_dpc1+p_dpc);
            p_fdsac = 1/2*(p_fdsac1+p_fdsac);
        end
        Result_ISAC(index,index_dB) = real(log2(det(eye(K)+pc/K*(H*H')))); % Equal Power Allocation
        Result_FDSAC(index,index_dB) = alpha*real(log2(det(eye(K)+pc/alpha/K*(H*H')))); % Equal Power Allocation
        Result_ISAC1(index,index_dB) = real(log2(det(eye(M)+(H'*diag(p_dpc)*H)))); % DPC-Based
        Result_FDSAC1(index,index_dB) = alpha*real(log2(det(eye(M)+1/alpha*(H'*diag(p_fdsac)*H)))); % DPC-Based
    end
    result(index) = real(log2(det(H*H')));
end
plot(CSNR_dB,mean(Result_ISAC),'-o');hold on;
plot(CSNR_dB,mean(Result_FDSAC),'-s');hold on;
plot(CSNR_dB,mean(Result_ISAC1),'--d');hold on;
plot(CSNR_dB,mean(Result_FDSAC1),'--s');hold on;
result = mean(result);
plot(CSNR_dB,K*log2(10.^(CSNR_dB/10)/K)+result);hold on; % High-SNR Slope
plot(CSNR_dB,alpha*(K*log2(10.^(CSNR_dB/10)/alpha/K)+result)); % High-SNR Slope
axis([-10,40,0,25]);
Data = [mean(Result_ISAC1);mean(Result_FDSAC1);K*log2(10.^(CSNR_dB/10)/K)+result;alpha*(K*log2(10.^(CSNR_dB/10)/alpha/K)+result);...
        mean(Result_ISAC);mean(Result_FDSAC)];
save Downlink_ISAC_ECR Data