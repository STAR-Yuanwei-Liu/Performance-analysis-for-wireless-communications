%% ----------------------------------------------------------------------------- %%
% This Matlab codes characterize the uplink outage probability of communication users in the paper:
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
Ps = 10; % Sensing SNR
s = waterfill(Ps,1./eig_t');
Radar_I = zeros(L,1);
Radar_I([1:1:M]) = s'.*eig_t; % Calculate the interference from sensing signals
Monte_Carlo = 1e6;
Outer_Monte_Carlo = 1e2;
CSNR_dB = [-10:5:40]; % Communication SNR
alpha = 0.5; % Spectrum allocation factor in FDSAC
R_target = 5; % Target communication rate
Result_ISAC = zeros(1,length(CSNR_dB));
Result_FDSAC = zeros(1,length(CSNR_dB));
L = 4; % Length of Radar Waveform
%% Calculate the communication rate
parfor outer_index = [1:1:Outer_Monte_Carlo]
    for index = [1:1:Monte_Carlo]
        [outer_index,index]
        H = 1/sqrt(2)*randn(K,N) + 1j*1/sqrt(2)*randn(K,N);
        Result_ISAC_Tmp = ones(1,length(CSNR_dB));
        Result_FDSAC_Tmp = ones(1,length(CSNR_dB));
        for index_dB = [1:1:length(CSNR_dB)]
            pc = 10^(CSNR_dB(index_dB)/10);
            Tmp = 0;
            for i = [1:1:L]
                sigma = 1 + Radar_I(i);
                Tmp = Tmp + (sign(R_target-real(log2(det(eye(K)+pc/sigma*(H*H')))))+1)/2;
            end
            Result_ISAC_Tmp(index_dB) = Tmp/L; % Communication rate in ISAC
            Result_FDSAC_Tmp(index_dB) = (sign(R_target-alpha*real(log2(det(eye(K)+pc/alpha*(H*H')))))+1)/2; % Communication rate in FDSAC
        end
        Result_ISAC = Result_ISAC + Result_ISAC_Tmp;
        Result_FDSAC = Result_FDSAC + Result_FDSAC_Tmp;
    end
end
semilogy(CSNR_dB,Result_ISAC/Outer_Monte_Carlo/Monte_Carlo,'-o');hold on;
semilogy(CSNR_dB,Result_FDSAC/Outer_Monte_Carlo/Monte_Carlo,'-s');hold on;
Data = [Result_ISAC/Outer_Monte_Carlo/Monte_Carlo;Result_FDSAC/Outer_Monte_Carlo/Monte_Carlo];
save Uplink_ISAC_OP Data