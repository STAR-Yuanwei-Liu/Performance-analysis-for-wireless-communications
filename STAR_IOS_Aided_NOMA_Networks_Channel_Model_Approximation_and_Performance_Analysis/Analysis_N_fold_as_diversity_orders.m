% The codes is for paper of "STAR-IOS Aided NOMA Networks: Channel Model Approximation and Performance Analysis"
%  in IEEE Transactions on Wireless Communications, vol. 21, no. 9, pp. 6861-6876, Sept. 2022, doi: 10.1109/TWC.2022.3152703.


clc
clear all;
%% setting
num = 1*1e4;          % simulation number             
P_dBM = 15:25;              %transmit power dBm
P_t= 10.^((P_dBM-30)./10);

BW=10*10^6; %10 MHz
Nf=10;%dB
sigma2_dbm= -170+10*log10(BW)+Nf; %Thermal noise in dBm
sigma_square=10^((sigma2_dbm-30)/10);

fc = 1e9; % 1 GHz carrier
c = 3*10^8; % speed of light
wavelength=3*10^8./fc;%wavelength
C_L=(wavelength/(4*pi))^2; %intercept of NLOS


rho = P_t./sigma_square; %%transmit power
rho_dB = 10.*log10(rho);

alpha = 2.4;                 % large-scale parameter

N=30; % num of RIS elements

k = 2; %% Rician parameter
s = sqrt(k./(1+k)); % matlab paramiter one
sigma = sqrt( 1./ (2.*(1+k))  ); % matlab paramiter Two

a_rfr = 0.6;
a_rfl = 0.4; %power allocation from BS

beta_rfr = 0.3; %power allocation by RIS
beta_rfl = 0.7;

R = 20; % RIS surving range
r_1 = 100; % distance between BS to RIS

R_user=0.1;  % threshold of SNR
R_SIC=0.1;
gamma_th =2.^(R_user)-1;
gamma_th_SIC=2.^(R_SIC)-1;


Gamma = max([gamma_th_SIC*sigma_square/(a_rfr-gamma_th_SIC*a_rfl), gamma_th*sigma_square/a_rfl]);
Gamma3 = gamma_th*sigma_square/(a_rfr-gamma_th*a_rfl);
k1 = k;
k2 = k;
%% diversity order
for i = 1:length(rho_dB)
miu = 4*sqrt(pi)*(1+k1)*(1+k2)/(exp(k1+k2)*gamma(2.5)) * hypergeom([1,1/2],5/2,1);

index1(i) = (Gamma3*r_1^alpha/(C_L^2*P_t(i)))^N;
P_out_rfr(i) = miu^N*R^(alpha*N)/(N*factorial(2*N-1)*(alpha*N+2)*beta_rfr^N) * index1(i) ;
index2(i) = (Gamma*r_1^alpha/(C_L^2*P_t(i)))^N;
P_out_rfl(i) = miu^N*R^(alpha*N)/(N*factorial(2*N-1)*(alpha*N+2)*beta_rfl^N) * index2(i) ;
end

%% order limit
P_inf_dBM = 110;              %transmit power dBm
P_inf= 10.^((P_inf_dBM-30)./10);
index1_inf = (Gamma3*r_1^alpha/(C_L^2*P_inf))^N;
P_out_rfr_inf = sqrt(beta_rfr)*miu^N*R^(alpha*N)/(N*factorial(2*N-1)*(alpha*N+2)) * index1_inf ;;
index2_inf = (Gamma*r_1^alpha/(C_L^2*P_inf))^N;
P_out_rfl_inf = sqrt(beta_rfl)*miu^N*R^(alpha*N)/(N*factorial(2*N-1)*(alpha*N+2)) * index2_inf ;

lim_rfl = -log(-P_out_rfl_inf)/log(P_inf);
lim_rfr = - log(-P_out_rfr_inf)/log(P_inf);
%% plot

semilogy(rho_dB,P_out_rfr,'ro-','LineWidth',2);
hold on;

semilogy(rho_dB,P_out_rfl,'b*-','LineWidth',2);
hold on;