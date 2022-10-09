% The codes for paper
% B. Zhao, C. Zhang, W. Yi and Y. Liu, "Ergodic Rate Analysis of STAR-RIS Aided NOMA Systems", IEEE Communications Letters; 2022.

clear
clc
%% Parameter setting
num = 1*1e4;
k=7;
theta=32.8;

P_dBM = 5:5:50; %transmit power in dBm
Pt= 10.^((P_dBM-30)./10);

a_far = 0.7;
a_near = 0.3; %power allocation from BS

H=50;
R1 = 400; % RIS surving range
dBR = 600; % distance between BS to RIS

alpha = 2.6;

BW=10*10^6; 
Nf=10;
sigma2_dbm= -170+10*log10(BW)+Nf; 
sigma_square=10^((sigma2_dbm-30)/10);

gamma_th_SIC=1;

rho=P_dBM-sigma2_dbm;

%% simulation
for i = 1:length(P_dBM)
    r_far = sqrt(R1^2.*rand(1,num));          
    r_near = sqrt(R1^2.*rand(1,num));
    index_far = zeros(1,num);
    index_near = zeros(1,num);

    for j =1:1:num        
        gamma_far(j)=random('Gamma',k,theta);
        gamma_near(j)=random('Gamma',k,theta);
        g_far(j) = gamma_far(j);
        g_near(j) =gamma_near(j);


        SNR_far(j) = a_far*Pt(i)*dBR^(-alpha)*(r_far(j)^2+H^2)^(-alpha/2)*g_far(j) / ...
                 (a_near*Pt(i)*dBR^(-alpha)*(r_far(j)^2+H^2)^(-alpha/2)*g_far(j) +sigma_square);

        SNR_SIC(j) = a_far*Pt(i)*dBR^(-alpha)*(r_near(j)^2+H^2)^(-alpha/2)*g_near(j) / ...
                 (a_near*Pt(i)*dBR^(-alpha)*(r_near(j)^2+H^2)^(-alpha/2)*g_near(j) +sigma_square);

        SNR_near(j) = a_near*Pt(i)*dBR^(-alpha)*(r_near(j)^2+H^2)^(-alpha/2)*g_near(j) / sigma_square;
             
         
         index_far(j) = log2(1+SNR_far(j));
                if SNR_SIC(j) > gamma_th_SIC 
                    index_near(j) = log2(1+SNR_near(j));
                else
                    index_near(j) = 0;
                end
    end
    
    Er_Far_ris(i) = sum(index_far)/num;
    Er_Near_ris(i) = sum(index_near)/num;
  
end

%% Plot
figure
plot(rho,Er_Far_ris,'r-');
ylim([0 3])
title('Far user')

figure
plot(rho,Er_Near_ris,'b-');
ylim([0 9])
title('Near user')
