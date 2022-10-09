clc;clear;close all;

%--------------------------------------------------------------------------
% "STAR-RIS Aided NOMA in Multicell Networks: A General Analytical Framework
% With Gamma Distributed Channel Modeling" in IEEE Transactions on communications
%
% Ziyi Xie
% Harbin Institute of Technology
%
% Main Code for Fig.2
%--------------------------------------------------------------------------

T = 2e6;
n = 751;  % Number of scale value

N = 4;    % number of elements in a STAR-RIS
K = 2;
sigma = sqrt(1/2);
m = 4;
omega = 1;
beta = 4;

%% simulation
% Rician channel
sim_rayl = raylrnd(sigma,T,N);
sim_NLoS = sum(sim_rayl,2);

sim_H1 = N*sqrt(K/(K+1)) + sqrt(1/(K+1))*sum(sim_rayl,2);
sim_H_p1 = sim_H1.^2;

[x1, sim_pdf1, sim_cdf1] = counting(sim_H_p1, min(sim_H_p1), 60, n);

% Rayleigh channel
sim_H2 = sum(sim_rayl,2);
sim_H_p2 = sim_H2.^2;

[x2, sim_pdf2, sim_cdf2] = counting(sim_H_p2, min(sim_H_p2), 45, n);

% Nakagami channel
pd3 = makedist('nakagami',m,omega);
sim_naka = sum(random(pd3,T,N),2);

sim_H3 = sum(sim_naka,2);
sim_H_p3 = sim_H3.^2;

[x3, sim_pdf3, sim_cdf3] = counting(sim_H_p3, min(sim_H_p3), 35, n);

% Weibull channel
pd4 = makedist('Weibull','a',omega,'b',beta);
sim_wei = sum(random(pd4,T,N),2);

sim_H4 = sum(sim_wei,2);
sim_H_p4 = sim_H4.^2;

[x4, sim_pdf4, sim_cdf4] = counting(sim_H_p4, min(sim_H_p4), 40, n);

% double-Rayleigh channel
sim_rayl2 = raylrnd(sigma,T,N);
sim_H5 = sum(sim_rayl.*sim_rayl2,2);
sim_H_p5 = sim_H5.^2;

[x5, sim_pdf5, sim_cdf5] = counting(sim_H_p5, min(sim_H_p5), 50, n);

% double-Rician Channel
% sim_rayl2 = raylrnd(sigma,T,N);
double_H = (sqrt(K/(K+1)) + sqrt(1/(K+1))*sim_rayl).*(sqrt(K/(K+1)) + sqrt(1/(K+1))*sim_rayl2);
sim_H6 = sum(double_H,2);
sim_H_p6 = sim_H6.^2;

[x6, sim_pdf6, sim_cdf6] = counting(sim_H_p6, min(sim_H_p6), 90, n);

%% Gamma Approximation (Analysis)
% Rician channel
m_ric = sigma*sqrt(pi/2)*sqrt(1/(K+1))+sqrt(K/(K+1));
v_ric = (4-pi)/2*sigma^2*(1/(K+1));

m_H_power1 = N^2*m_ric^2 + N*v_ric;
v_H_power1 = 4*m_ric^2*v_ric*N^3 + 2*v_ric^2*N^2;

shape1 = m_H_power1^2 / v_H_power1;   % shapex/scalex:  eq.(12)
scale1 = v_H_power1 / m_H_power1;

shape11 = m_ric^2/v_ric/4*N;          % shapexx/scalexx:  eq.(13)  
scale11 = 4*v_ric*N;

z = 1:50:n;
approx_cdf1 = gammainc(x1(z)/scale1,shape1);
approx_cdf11 = gammainc(x1(z)/scale11,shape11);

% Rayleigh channel
m_rayl = sigma*sqrt(pi/2);
v_rayl = (4-pi)/2*sigma^2;

m_H_power2 = N^2*m_rayl^2 + N*v_rayl;
v_H_power2 = 4*m_rayl^2*v_rayl*N^3 + 2*v_rayl^2*N^2;

shape2 = m_H_power2^2 / v_H_power2;
scale2 = v_H_power2 / m_H_power2;

shape22 = m_rayl^2/v_rayl/4*N;
scale22 = 4*v_rayl*N;

approx_cdf2 = gammainc(x2(z)/scale2,shape2);
approx_cdf22 = gammainc(x2(z)/scale22,shape22);

% Nakagami channel
m_naka = gamma(m+0.5)/gamma(m)*(omega/m)^0.5;
v_naka = omega-omega/m*(gamma(m+0.5)/gamma(m))^2;

m_H_power3 = N^2*m_naka^2 + N*v_naka;
v_H_power3 = 4*m_naka^2*v_naka*N^3 + 2*v_naka^2*N^2;

shape3 = m_H_power3^2 / v_H_power3;
scale3 = v_H_power3 / m_H_power3;

shape33 = m_naka^2/v_naka/4*N;
scale33 = 4*v_naka*N;

approx_cdf3 = gammainc(x3(z)/scale3,shape3);
approx_cdf33 = gammainc(x3(z)/scale33,shape33);

% Weibull channel
m_wei = omega*gamma(1+1/beta);
v_wei = omega^2*(gamma(1+2/beta)-(gamma(1+1/beta))^2);

m_H_power4 = N^2*m_wei^2 + N*v_wei;
v_H_power4 = 4*m_wei^2*v_wei*N^3 + 2*v_wei^2*N^2;

shape4 = m_H_power4^2 / v_H_power4;
scale4 = v_H_power4 / m_H_power4;

shape44 = m_wei^2/v_wei/4*N;
scale44 = 4*v_wei*N;

approx_cdf4 = gammainc(x4(z)/scale4,shape4);
approx_cdf44 = gammainc(x4(z)/scale44,shape44);

% double-Rayleigh channel
m_drayl = sigma^2*pi/2;
v_drayl = 4*(1-pi^2/16)*sigma^4;

m_H_power5 = N^2*m_drayl^2 + N*v_drayl;
v_H_power5 = 4*m_drayl^2*v_drayl*N^3 + 2*v_drayl^2*N^2;

shape5 = m_H_power5^2 / v_H_power5;
scale5 = v_H_power5 / m_H_power5;

shape55 = m_drayl^2/v_drayl/4*N;
scale55 = 4*v_drayl*N;

approx_cdf5 = gammainc(x5(z)/scale5,shape5);
approx_cdf55 = gammainc(x5(z)/scale55,shape55);

% double-Rician Channel
a = sqrt(1/(K+1));c = sqrt(K/(K+1));
xx2 = 2*sigma^2;xx1 = sqrt(pi/2)*sigma;
m_dric = a^2*xx1^2 + 2*a*c*xx1 + c^2;
v_dric = a^4*xx2^2 + 4*a^3*c*xx1*xx2 + 2*a^2*c^2*xx1^2 + 2*a^2*c^2*(xx2+xx1^2) + 4*a*c^3*xx1 +c^4 - m_dric^2;

m_H_power6 = N^2*m_dric^2 + N*v_dric;
v_H_power6 = 4*m_dric^2*v_dric*N^3 + 2*v_dric^2*N^2;

shape6 = m_H_power6^2 / v_H_power6;
scale6 = v_H_power6 / m_H_power6;

shape66 = m_dric^2/v_dric/4*N;
scale66 = 4*v_dric*N;

approx_cdf6 = gammainc(x6(z)/scale6,shape6);
approx_cdf66 = gammainc(x6(z)/scale66,shape66);

%% plot
% Rician channel
s1 = plot(x1,sim_cdf1,'-','color',[1,0.6,0.6],'linewidth',2);hold on;
a1 = plot(x1(z),approx_cdf1,'*','MarkerSize',6,'color',[0,0,0],'linewidth',0.6);hold on;
a11 = plot(x1(z),approx_cdf11,'o','color',[0.3,0.3,0.3],'linewidth',0.6);hold on;

% Rayleigh channel
s2 = plot(x2,sim_cdf2,'-','color',[0,0,0],'linewidth',1);hold on;
plot(x2(z),approx_cdf2,'*','MarkerSize',6,'color',[0,0,0],'linewidth',0.6);hold on;
plot(x2(z),approx_cdf22,'o','color',[0.3,0.3,0.3],'linewidth',0.6);hold on;

% Nakagami channel
s3 = plot(x3,sim_cdf3,'-.','color',[1,0,0],'linewidth',1);hold on;
plot(x3(z),approx_cdf3,'*','MarkerSize',6,'color',[0,0,0],'linewidth',0.6);hold on;
plot(x3(z),approx_cdf33,'o','color',[0.3,0.3,0.3],'linewidth',0.6);hold on;

% Weibull channel
s4 = plot(x4,sim_cdf4,':','color',[0,0,1],'linewidth',1.1);hold on;
plot(x4(z),approx_cdf4,'*','MarkerSize',6,'color',[0,0,0],'linewidth',0.6);hold on;
plot(x4(z),approx_cdf44,'o','color',[0.3,0.3,0.3],'linewidth',0.6);hold on;

% double-Rayleigh channel
s5 = plot(x5,sim_cdf5,'--','color',[0,0,0],'linewidth',1);hold on;
plot(x5(z),approx_cdf5,'*','MarkerSize',6,'color',[0,0,0],'linewidth',0.6);hold on;
plot(x5(z),approx_cdf55,'o','color',[0.3,0.3,0.3],'linewidth',0.6);hold on;

% double-Rayleigh channel
s6 = plot(x6,sim_cdf6,'-.','color',[0.1,0.9,0.1],'linewidth',2);hold on;
plot(x6(z),approx_cdf6,'*','MarkerSize',6,'color',[0,0,0],'linewidth',0.6);hold on;
plot(x6(z),approx_cdf66,'o','color',[0.3,0.3,0.3],'linewidth',0.6);hold on;

legend([s2 s3 s1 s4 s5 s6 a1 a11],'Rayleigh Channel','Nakagami Channel','Rician Channel','Weibull Channel',...
    'Double-Rayleigh Channel','Double-Rician Channel','Lemma 1','Corollary 1','FontName','Times New Roman','FontSize',10);

xlabel('|{\it h_{r,S}}|^2','FontName','Times New Roman');
ylabel('CDF of |{\it h_{r,S}}|^2','FontName','Times New Roman');