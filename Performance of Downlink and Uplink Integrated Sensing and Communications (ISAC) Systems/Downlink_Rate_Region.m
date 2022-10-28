%% ----------------------------------------------------------------------------- %%
% This Matlab codes characterize the downlink rate region in the paper:
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
Ps = 10;
resolution = [0:0.001:1];
FDSAC_ECR = zeros(1,length(resolution));
FDSAC_SR = zeros(1,length(resolution));
ISAC_ECR = zeros(1,length(resolution));
ISAC_SR = zeros(1,length(resolution));
Monte_Carlo = 1e4;
%% Calculate the sensing rate of FDSAC
for index = [1:1:length(resolution)]
    alpha = resolution(index);
    ps = 10^(Ps/10);
    if alpha~=1
        s = waterfill(ps,(1-alpha)./(eig_t'));
        FDSAC_SR(index) = (1-alpha)*N/L*sum(log2(1+1/(1-alpha)*s.*(eig_t')));
    end
end
%% Calculate the communication rate in ISAC
%% Calculate the communiction rate in FDSAC
Interference = zeros(1,length(resolution));
parfor C_index = [1:1:Monte_Carlo]
    H = (1/sqrt(2)*randn(K,M) + 1j*1/sqrt(2)*randn(K,M))*(R^0.5); % Correlated Rayleigh Fading
    h1 = H(1,:); % Channel of communication user 1
    h2 = H(2,:); % Channel of communication user 2
    v = real([h1*h1',h2*h2']);
    iteration = 20;
    %% Iterative Water-Filling Method
    % For more details, please refer to the following paper
    % Nihar Jindal, Wonjong Rhee, Sriram Vishwanath, Syed Ali Jafar, and Andrea Goldsmith, 
    % "Sum Power Iterative Water-filling for Multi-Antenna Gaussian Broadcast Channels," 
    % IEEE Transactions on Information Theory, Vol. 51, No. 4, pp. 1570-1580, April 2005.
    ISAC_ECR_Tmp = zeros(1,length(resolution));
    FDSAC_ECR_Tmp = zeros(1,length(resolution));
    Interference_Tmp = zeros(1,length(resolution));
    for index = [1:1:length(resolution)]
        [index,C_index]
        alpha = resolution(index);
        if alpha~=0
            pc = alpha*(10^(Pc/10));
            p_dpc = pc/K*[1,1];
            p_fdsac = pc/K*[1,1];
            for j = [1:1:iteration]    
                p_dpc1 = waterfill(pc,1./(v - (p_dpc*((abs(h1*h2'))^2))./(1+p_dpc.*wrev(v))));
                p_dpc = 1/2*(p_dpc1+p_dpc);
                p_fdsac1 = waterfill(pc,1./(v - (p_fdsac*((abs(h1*h2'))^2))./(1+p_fdsac.*wrev(v)))*alpha);
                p_fdsac = 1/2*(p_fdsac1+p_fdsac);
            end
            ISAC_ECR_Tmp(index) = real(log2(det(eye(M)+(H'*diag(p_dpc)*H)))); % Communication rate in ISAC
            FDSAC_ECR_Tmp(index) = alpha*real(log2(det(eye(M)+1/alpha*(H'*diag(p_fdsac)*H)))); % Communication rate in FDSAC
            B1 = eye(M) - h2'*h2*p_dpc(2)/(1+p_dpc(2)*h2*h2');
            M1 = B1*(h1'*h1)*B1*p_dpc(1)/(h1*B1*h1');
            M2 = (h2'*h2)*p_dpc(2)*(1+h2*M1*h2')/(h2*h2');
            Interference_Tmp(index) = real(trace(Rt*(M1 + M2)));
        end
    end
    ISAC_ECR = ISAC_ECR + ISAC_ECR_Tmp;
    FDSAC_ECR = FDSAC_ECR + FDSAC_ECR_Tmp;
    Interference = Interference + Interference_Tmp;
end
Interference = Interference/Monte_Carlo;
%% Calculate the sensing rate in ISAC
for index = [1:1:length(resolution)]
    tmp = Interference(index);
    s = waterfill(ps,(1 + tmp)./eig_t');
    ISAC_SR(index) = N/L*sum(log2(1+1/(1 + tmp)*s.*(eig_t')));
end
FDSAC_ECR = FDSAC_ECR/Monte_Carlo; 
ISAC_ECR = ISAC_ECR/Monte_Carlo;
plot(ISAC_SR,ISAC_ECR,':');
hold on;
plot(FDSAC_SR,FDSAC_ECR,'-.k');
hold on;
plot([0,ISAC_SR(end)],[ISAC_ECR(end),ISAC_ECR(end)],'k-');
hold on;
plot([ISAC_SR(end),ISAC_SR(end)],[0,ISAC_ECR(end)],'k-');
hold on;
Data = [ISAC_SR;ISAC_ECR;FDSAC_SR;FDSAC_ECR];
save Downlink_Rate_Region Data