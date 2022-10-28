%% ----------------------------------------------------------------------------- %%
% This Matlab codes characterize the rate region in the paper:
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
% kappa = 0.5;% Power allocation factor in FDSAC
% mu = 0.5;% Spectrum allocation factor in FDSAC
eig = [1;0.1;0.05;0.01];% Eigenvalues of the correlation matrix of the target response matrix
R = diag(eig);
Monte = 1e2;
CR_SR = cell(1,5);
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
SNR_dB = [5];%Sensing SNR
ps = 10^(SNR_dB/10);
%% Sensing-Centric
s = waterfill(ps,1./eig'/L);
CR_SC = mean(sum(transpose(log2(1+(ones(Monte,1)*s).*H))));
SR_SC = N/L*sum(log2(1+s.*(eig')*L));
CR = [0,CR_SC,CR_SC]; % Calculate the Communication Rate
SR = [SR_SC,SR_SC,0]; % Calculate the Sensing Rate
plot(CR,SR,'-.'); hold on;
CR_SR{1,1}=[CR_SC,SR_SC];
%% Communication-Centric
s = waterfill(ps*ones(Monte,1),1./H);
CR_CC = mean(sum(transpose(log2(1+s.*H))));
SR_CC = mean(real(N/L*sum(transpose(log2(1+s.*(ones(Monte,1)*(eig'))*L)))));
CR = [0,CR_CC,CR_CC]; % Calculate the Communication Rate
SR = [SR_CC,SR_CC,0]; % Calculate the Sensing Rate
plot(CR,SR,'--'); hold on;
CR_SR{1,2}=[CR_CC,SR_CC];
%% Equal Power Allocation
s = ps/M*ones(1,M);
CR_EPA = mean(sum(transpose(log2(1+(ones(Monte,1)*s).*H))));
SR_EPA = N/L*sum(log2(1+s.*(eig')*L)); 
CR = [0,CR_EPA,CR_EPA]; % Calculate the Communication Rate
SR = [SR_EPA,SR_EPA,0]; % Calculate the Sensing Rate
plot(CR,SR,'--'); hold on;
CR_SR{1,3}=[CR_EPA,SR_EPA];
%% FDSAC
precision = 1e-2;
resolution = [0:precision:1];
FDSAC_Rate_Region = ones(2,length(resolution)*length(resolution));
for a = [1:1:length(resolution)]
    for b = [1:1:length(resolution)]
        [a,b]
        mu = resolution(a); % Power Allocation Factor
        kappa = resolution(b); % Spectrum Allocation Factor
        if kappa == 1
            SR = 0;
        else
            s = waterfill(ps*(1-mu),1./eig'/L/(1-kappa));
            SR = N/L*(1-kappa)*sum(log2(1+s.*(eig')*L/(1-kappa))); % Calculate the Sensing Rate          
        end
        if kappa == 0
            CR = 0;
        else
            s = waterfill(ps*mu,1./H/kappa);
            CR = kappa*mean(sum(transpose(log2(1+s.*H/kappa)))); % Calculate the Communication Rate
        end
        FDSAC_Rate_Region(1,b+(a-1)*length(resolution)) = SR;
        FDSAC_Rate_Region(2,b+(a-1)*length(resolution)) = CR;
    end
end
x = (FDSAC_Rate_Region(2,:))';
y = (FDSAC_Rate_Region(1,:))';
dt = delaunayTriangulation(x,y);
k = convexHull(dt); % Calculate the Convex Hull of all Communication-Sensing Rate Pairs
plot(dt.Points(k,1),dt.Points(k,2), '-');hold on;
CR_SR{1,4}=[dt.Points(k,1);dt.Points(k,2)];
%% Trade-Off
eig = [1,0.1,0.05,0.01];
precision = 1e-2; % Precision of Rate Profile
rate_profile_factor = [precision:precision:(1-precision)]; % Rate Profile Factor
CR = zeros(1,length(rate_profile_factor)); % Used to Storage the Communication Rate
SR = zeros(1,length(rate_profile_factor)); % Used to Storage the Sensing Rate
for outer_index = [1:1:Monte]
    CR_tmp = ones(1,length(rate_profile_factor)); % Used to Storage the Communication Rate
    SR_tmp = ones(1,length(rate_profile_factor)); % Used to Storage the Sensing Rate
    for index = [1:1:length(rate_profile_factor)]
        index
        varepsilon = rate_profile_factor(index); % Rate Profile Factor
        cvx_clear % CVX is used to solve the rate profile problem
        cvx_begin quiet
        % State the Variables. Note we must highlight the size and the semidefinition of W £¨complex semidefinite£©
        variable s(1,M) % The power allocation
        variable R  % The second variable: Rate R
        minimize(-1*real(R)); % Objective function
        % Note that we must use the function to extract the real part in CVX programming
        % The above rule also applies to the programming of the constraints
        subject to % The follow statements are the constraints of the original problem
        (R*(1-varepsilon)) - sum_log(1+s.*eig*L)/log(2) <= 0 % Rate Constraint of Sensing
        (R*varepsilon) - sum_log(1+s.*H(outer_index,:))/log(2) <= 0; % Rate Constraint of Communication
        real(sum(s)) - ps <= 0; % Power Constraint
        0 <= s; % Power Constraint
        cvx_end
        CR_tmp(index) = sum(log2(1+s.*H(outer_index,:))); % Calculate the Communication Rate
        SR_tmp(index) = real(N/L*sum(log2(1+s.*(eig)*L))); % Calculate the Sensing Rate
    end
    CR = CR + CR_tmp;
    SR = SR + SR_tmp;
end
CR_Pareto = CR/Monte;
SR_Pareto = SR/Monte;
plot(CR_Pareto,SR_Pareto);hold on;
CR_SR{1,5}=[CR_Pareto;SR_Pareto];
save Rate_Region_Data CR_SR