% The codes is for paper of "STAR-IOS Aided NOMA Networks: Channel Model Approximation and Performance Analysis"
%  in IEEE Transactions on Wireless Communications, vol. 21, no. 9, pp. 6861-6876, Sept. 2022, doi: 10.1109/TWC.2022.3152703.

clc
clear all;
%% setting
num = 1*1e5;          % simulation number             
P_dBM = 15:2:25;              %transmit power dBm
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

%% simulation
P_out_sim_rfr = zeros(1,length(rho));
P_out_sim_rfl = zeros(1,length(rho));
for i = 1:length(rho_dB)
    tic;
    r_rfr = sqrt(R^2.*rand(1,num));          % user position :reflection
    r_rfl = sqrt(R^2.*rand(1,num));          % user position :refraction
    Poutsum_rfr = 0;
    Poutsum_rfl = 0;
    for j =1:1:num
    rician1 = random('rician',s,sigma,[1,N]);
    rician2 = random('rician',s,sigma,[1,N]);
    g_rfr = beta_rfr.*sum(rician1.*rician2).^2;
    g_rfl = beta_rfl.*sum(rician1.*rician2).^2;
    
    SNR_SIC(j) = a_rfr*P_t(i)*C_L*r_1^(-alpha)*C_L*r_rfl(j)^(-alpha)*g_rfl / ...
                 (a_rfl*P_t(i)*C_L*r_1^(-alpha)*C_L*r_rfl(j)^(-alpha)*g_rfl +sigma_square);
    SNR_rfl(j) = a_rfl*P_t(i)*C_L*r_1^(-alpha)*C_L*r_rfl(j)^(-alpha)*g_rfl /sigma_square ;
    SNR_rfr(j) = a_rfr*P_t(i)*C_L*r_1^(-alpha)*C_L*r_rfr(j)^(-alpha)*g_rfr / ...
                 (a_rfl*P_t(i)*C_L*r_1^(-alpha)*C_L*r_rfr(j)^(-alpha)*g_rfr +sigma_square);
             
         if   SNR_SIC(j) < gamma_th_SIC   || SNR_rfl(j) < gamma_th 
             Poutsum_rfl= Poutsum_rfl+1;
         end
         if SNR_rfr(j) < gamma_th
             Poutsum_rfr = Poutsum_rfr +1;
         end
     
    end
     P_out_sim_rfr(i) = Poutsum_rfr/num;
     P_out_sim_rfl(i) = Poutsum_rfl/num;
    
    toc;
end

% figure
semilogy(rho_dB,P_out_sim_rfr,'ro-','LineWidth',2);
hold on;

semilogy(rho_dB,P_out_sim_rfl,'b*-','LineWidth',2);
hold on;
