% The codes for paper
% B. Zhao, C. Zhang, W. Yi and Y. Liu, "Ergodic Rate Analysis of STAR-RIS Aided NOMA Systems", IEEE Communications Letters; 2022.

clear
clc
%% Parameter setting
num =1*1e4;
k=3;
theta=14;
N=30;
           
P_dBM = 5:5:50; %transmit power dBm 0-40
Pt= 10.^((P_dBM-30)./10);

beta_far = 0.5; 
beta_near = 0.5; %power allocation by STAR-RIS

a_far = 0.7;
a_near = 0.3; %power allocation from BS
a_index = a_far/a_near;

H=30;
R1 = 100;
R2 = 200;    
dBR = 400; % distance between BS to STAR-RIS

alpha = 2.6;

BW=10*10^6; 
Nf=10;
sigma2_dbm= -170+10*log10(BW)+Nf; %Thermal noise in dBm
sigma_square=10^((sigma2_dbm-30)/10);

gamma_th_SIC=1;

rho=P_dBM-sigma2_dbm;

%% simulation

% 5x6
x=[1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5];
y=[1,2,3,4,5,6];
y=[y,y,y,y,y];

% 5x10
% x=[1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5];
% y=[1,2,3,4,5,6,7,8,9,10];
% y=[y,y,y,y,y];

% 7x10
% x=[1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7];
% y=[1,2,3,4,5,6,7,8,9,10];
% y=[y,y,y,y,y,y,y];

% correlation matrix
R=zeros(N);
for p=1:N
    for q=1:N
        a=[x(p),y(p)];
        b=[x(q),y(q)];
        R(p,q)=0.75*sinc(norm(a-b));
    end
end


for i = 1:length(P_dBM)
    tic;
    r_far = sqrt(R1^2+(R2^2-R1^2).*rand(1,num));       
    r_near = sqrt(R1^2.*rand(1,num));
    index_far = zeros(1,num);
    index_near = zeros(1,num);
    index_far_co = zeros(1,num);
    index_near_co = zeros(1,num);
    for j =1:1:num
        gamma_far(j)=random('Gamma',k,theta);
        gamma_near(j)=random('Gamma',k,theta);
        g_far(j) = beta_far.*gamma_far(j);
        g_near(j) = beta_near.*gamma_near(j);

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

        g_far_co(j) = trace(R^0.5*(beta_far.*eye(N))*R*(beta_far.*eye(N))*R^0.5);
        g_near_co(j) = trace(R^0.5*(beta_near.*eye(N))*R*(beta_near.*eye(N))*R^0.5);

        SNR_far_co(j) = a_far*Pt(i)*dBR^(-alpha)*(r_far(j)^2+H^2)^(-alpha/2)*g_far_co(j) / ...
                 (a_near*Pt(i)*dBR^(-alpha)*(r_far(j)^2+H^2)^(-alpha/2)*g_far_co(j) +sigma_square);

        SNR_SIC_co(j) = a_far*Pt(i)*dBR^(-alpha)*(r_near(j)^2+H^2)^(-alpha/2)*g_near_co(j) / ...
                 (a_near*Pt(i)*dBR^(-alpha)*(r_near(j)^2+H^2)^(-alpha/2)*g_near_co(j) +sigma_square);

        SNR_near_co(j) = a_near*Pt(i)*dBR^(-alpha)*(r_near(j)^2+H^2)^(-alpha/2)*g_near_co(j) / sigma_square;
             
         
                index_far_co(j) = log2(1+SNR_far_co(j));
                if SNR_SIC_co(j) > gamma_th_SIC 
                    index_near_co(j) = log2(1+SNR_near_co(j));
                else
                    index_near_co(j) = 0;
                end

    end
    
    ErFar(i) = sum(index_far)/num;
    ErNear(i) = sum(index_near)/num;
    ErFar_co(i) = sum(index_far_co)/num;
    ErNear_co(i) = sum(index_near_co)/num;
    toc;
end

%% plot
figure
plot(rho,ErFar,'ro-' );
hold on;
plot(rho,ErFar_co,'k*-' );
hold on;
legend('uncorrelated','correlated')
ylim([0 3])
title('Far user')

figure
plot(rho,ErNear,'bo-' );
hold on;
plot(rho,ErNear_co,'k*-');
hold on;
legend('uncorrelated','correlated')
ylim([0 9])
title('Near user')
