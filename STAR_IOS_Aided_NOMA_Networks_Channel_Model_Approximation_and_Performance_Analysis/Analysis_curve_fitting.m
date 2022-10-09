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
% 
A= 30;
B_rfr = 22.46;% 6.739;
B_rfl = 22.46;% 15.58 ;


% A= 20;
% B_rfr = 15.06;
% B_rfl = B_rfr;

Gamma = max([gamma_th_SIC*sigma_square/(a_rfr-gamma_th_SIC*a_rfl), gamma_th*sigma_square/a_rfl]);
Gamma3 = gamma_th*sigma_square/(a_rfr-gamma_th*a_rfl);


%% analysis with curve fitting
for i = 1:length(rho)
    if a_rfr >gamma_th_SIC*a_rfl && a_rfr >gamma_th*a_rfl
        index_rfr = Gamma3/(P_t(i)*C_L^2)*r_1^alpha/B_rfr/beta_rfr ;
        index_rfl = Gamma/(P_t(i)*C_L^2)*r_1^alpha/B_rfl/beta_rfl ;
        P_out_rfl(i) = 2/R^2*integral(@(x)x.*gammainc(index_rfl.*x.^alpha,A) ,0,R);
        P_out_rfr(i) = 2/R^2*integral(@(x)x.*gammainc(index_rfr.*x.^alpha,A) ,0,R);
    end
end

semilogy(rho_dB,P_out_rfl,'r-')
hold on;
semilogy(rho_dB,P_out_rfr,'b-')
hold on;
