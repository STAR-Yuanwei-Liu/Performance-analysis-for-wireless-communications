%% ----------------------------------------------------------------------------- %%
% This Matlab codes characterize the uplink ergodic rate of communication users in the paper:
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
Monte_Carlo = 5e4;
CSNR_dB = [-10:5:50]; % Communication SNR
Result_ISAC = ones(Monte_Carlo,length(CSNR_dB));
Result_FDSAC = ones(Monte_Carlo,length(CSNR_dB));
alpha = 0.5; % Spectrum allocation factor in FDSAC
%% Calculate the communication rate
for index = [1:1:Monte_Carlo]
    index
    H = 1/sqrt(2)*randn(K,N) + 1j*1/sqrt(2)*randn(K,N);
    for index_dB = [1:1:length(CSNR_dB)]
        pc = 10^(CSNR_dB(index_dB)/10);
        Tmp = 0;
        for i = [1:1:L]
            sigma = 1 + Radar_I(i);
            Tmp = Tmp + real(log2(det(eye(K)+pc/sigma*(H*H'))));
        end
        Result_ISAC(index,index_dB) = 1/L*Tmp; % Communication rate in ISAC
        Result_FDSAC(index,index_dB) = alpha*real(log2(det(eye(K)+pc/alpha*(H*H')))); % Communication rate in FDSAC
    end
end
plot(CSNR_dB,mean(Result_ISAC),'-o');hold on;
plot(CSNR_dB,mean(Result_FDSAC),'-s');hold on;
%% Calculate the asymptotic communication rate
result = -K*double(eulergamma);
for l = [0:1:(K-1)]
    for i = [1:1:(N-l-1)]
        result = result + 1/i;
    end
end
Asym_ISAC = 0;
for i = [1:1:L]
    sigma = 1 + Radar_I(i);
    Asym_ISAC = Asym_ISAC + K*log2(10.^(CSNR_dB/10)/sigma)+result/log(2);
end
plot(CSNR_dB,Asym_ISAC/L);hold on;
plot(CSNR_dB,alpha*(K*log2(10.^(CSNR_dB/10)/alpha)+result/log(2)));
axis([-10,50,0,25]);
Data = [mean(Result_ISAC);mean(Result_FDSAC);Asym_ISAC/L;alpha*(K*log2(10.^(CSNR_dB/10)/alpha)+result/log(2))];
save Uplink_ISAC_ECR Data