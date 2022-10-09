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

% center limit parameters
E_ana = sqrt(pi/4/(1+k))*hypergeom(-0.5,1,-k);
V_ana = 1-(pi/4/(1+k))* (hypergeom(-0.5,1,-k))^2;

Heq_rfr = sqrt(beta_rfr)*N* E_ana^2;
Veq_rfr = beta_rfr*N*(2*E_ana^2*V_ana+V_ana^2);

Heq_rfl = sqrt(beta_rfl)*N* E_ana^2;
Veq_rfl = beta_rfl*N*(2*E_ana^2*V_ana+V_ana^2);

Gamma = max([gamma_th_SIC*sigma_square/(a_rfr-gamma_th_SIC*a_rfl), gamma_th*sigma_square/a_rfl]);
Gamma3 = gamma_th*sigma_square/(a_rfr-gamma_th*a_rfl);

%% Analytical results for center limit model 
for i = 1:length(rho)
    if a_rfr >gamma_th_SIC*a_rfl
        index = ( Gamma/(P_t(i)*C_L^2) )^0.5;
        fun = @(x) x.* (      erf( (Heq_rfl+index.*x.^(alpha/2).*r_1.^(alpha/2))/(sqrt(2*Veq_rfl) ) ) -...
           erf( (Heq_rfl-index.*x.^(alpha/2).*r_1.^(alpha/2))/(sqrt(2*Veq_rfl) ) )        );   
        P_out_rfl(i) = 1/R^2*integral(fun,0,R);
        
        index2 = ( Gamma3/(P_t(i)*C_L^2) )^0.5;
        fun2 = @(x) x.* (      erf( (Heq_rfr+index2.*x.^(alpha/2).*r_1.^(alpha/2))/(sqrt(2*Veq_rfr) ) ) -...
           erf( (Heq_rfr-index2.*x.^(alpha/2).*r_1.^(alpha/2))/(sqrt(2*Veq_rfr) ) ));
       P_out_rfr(i) = 1/R^2*integral(fun2,0,R);
    end
end
semilogy(rho_dB,P_out_rfl,'b--')
hold on;
semilogy(rho_dB,P_out_rfr,'r--')
hold on;
