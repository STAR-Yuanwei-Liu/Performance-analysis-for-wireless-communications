%% ----------------------------------------------------------------------------- %%
% This Matlab codes characterize the downlink outage probability of communication users in the paper:
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
Monte_Carlo = 1e10; % This should be set to a very large value to make the curve smooth
Outer_Monte_Carlo = 1e4; % This should be set to a very large value to make the curve smooth
CSNR_dB = [-10:5:40]; % SNR
Result_ISAC = zeros(1,length(CSNR_dB));
Result_FDSAC = zeros(1,length(CSNR_dB));
Result_ISAC1 = zeros(1,length(CSNR_dB));
Result_FDSAC1 = zeros(1,length(CSNR_dB));
alpha = 0.5; % Spectrum allocation factor in FDSAC
R_target = 5; % Target rate of sum communication rate
parfor outer_index = [1:1:Outer_Monte_Carlo]
    for index = [1:1:Monte_Carlo]
        [outer_index,index]
        H = (1/sqrt(2)*randn(K,M) + 1j*1/sqrt(2)*randn(K,M))*(R^0.5); % Correlated Rayleigh Fading
        h1 = H(1,:); % Channel of communication user 1
        h2 = H(2,:); % Channel of communication user 2
        v = real([h1*h1',h2*h2']);
        %% Iterative Water-Filling Method
        % For more details, please refer to the following paper
        % Nihar Jindal, Wonjong Rhee, Sriram Vishwanath, Syed Ali Jafar, and Andrea Goldsmith, 
        % "Sum Power Iterative Water-filling for Multi-Antenna Gaussian Broadcast Channels," 
        % IEEE Transactions on Information Theory, Vol. 51, No. 4, pp. 1570-1580, April 2005.
        iteration = 40;
        Result_ISAC_Tmp = ones(1,length(CSNR_dB));
        Result_FDSAC_Tmp = ones(1,length(CSNR_dB));
        Result_ISAC1_Tmp = ones(1,length(CSNR_dB));
        Result_FDSAC1_Tmp = ones(1,length(CSNR_dB));
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
            Result_ISAC_Tmp(index_dB) = (sign(R_target-real(log2(det(eye(K)+pc/K*(H*H')))))+1)/2; % Equal Power Allocation
            Result_FDSAC_Tmp(index_dB) = (sign(R_target-alpha*real(log2(det(eye(K)+pc/alpha/K*(H*H')))))+1)/2; % Equal Power Allocation
            Result_ISAC1_Tmp(index_dB) = (sign(R_target-real(log2(det(eye(M)+(H'*diag(p_dpc)*H)))))+1)/2; % DPC-Based
            Result_FDSAC1_Tmp(index_dB) = (sign(R_target-alpha*real(log2(det(eye(M)+1/alpha*(H'*diag(p_fdsac)*H)))))+1)/2; % DPC-Based
        end
        Result_ISAC = Result_ISAC + Result_ISAC_Tmp;
        Result_FDSAC = Result_FDSAC + Result_FDSAC_Tmp;
        Result_ISAC1 = Result_ISAC1 + Result_ISAC1_Tmp;
        Result_FDSAC1 = Result_FDSAC1 + Result_FDSAC1_Tmp;
    end
end
semilogy(CSNR_dB,Result_ISAC/Monte_Carlo/Outer_Monte_Carlo,'-o');hold on;
semilogy(CSNR_dB,Result_FDSAC/Monte_Carlo/Outer_Monte_Carlo,'-s');hold on;
semilogy(CSNR_dB,Result_ISAC1/Monte_Carlo/Outer_Monte_Carlo,'--d');hold on;
semilogy(CSNR_dB,Result_FDSAC1/Monte_Carlo/Outer_Monte_Carlo,'--s');hold on;
Data = [Result_ISAC;Result_FDSAC;Result_ISAC1;Result_FDSAC1]/Monte_Carlo/Outer_Monte_Carlo;
save Downlink_ISAC_OP Data