% The codes for paper
% B. Zhao, C. Zhang, W. Yi and Y. Liu, "Ergodic Rate Analysis of STAR-RIS Aided NOMA Systems", IEEE Communications Letters; 2022.

clear
clc
%% Parameter setting
num =1*1e4;
k=3;
theta=14;
           
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


%% analytical

for ii = 1:length(P_dBM)
    Phy = @(y) Pt(ii).*(y.^2+H^2).^(-0.5*alpha).*dBR^(-alpha);

    %% Far user
    t=@(x) sigma_square*x./(a_near*beta_near*Pt(ii)*dBR^(-alpha)*theta*(a_index-x));
    for jj = 1:k
        n = jj-1;
        
        fun=@(x) t(x).^(-2/alpha)./(x+1).*(gammainc(t(x)*(R2^2+H^2)^(alpha/2),n+2/alpha)*gamma(n+2/alpha)-gammainc(t(x)*(R1^2+H^2)^(alpha/2),n+2/alpha)*gamma(n+2/alpha));
     
       int_far(jj) = 1/factorial(n)*integral(@(x) fun(x),0,a_index);
      
    end
     R_far(ii) = 2/(log(2)*alpha*(R2^2-R1^2))*sum(int_far);

     %% Near user
     a_cal = a_near*beta_near*gamma_th_SIC/(a_far*beta_far-gamma_th_SIC*a_near*beta_near);
     ee=@(y) sigma_square./(a_near*beta_near.*Phy(y).*theta);
     aa=@(y) ee(y).*(a_cal+1);
     func_ori = @(y) log(a_cal+1)/log(2)*(1-gammainc(sigma_square.*gamma_th_SIC./(a_near.*beta_near.*Phy(y).*theta.*(a_index-gamma_th_SIC)),k, 'lower')).*2.*y./R1^2;
   
     fun1= @(y) exp(ee(y)).*expint(aa(y)).*y.*2./(R1^2);
       
    for jj=1:k-1
        n=jj;
        int=zeros(1,n);
        f2=@(x,y) (ee(y).*x).^n./(1+x).*exp(-ee(y).*x).*y.*2./(R1^2);
        fun31= @(y) ee(y).^n.*exp(ee(y)).*(-1)^n.*expint(aa(y)).*y.*2./(R1^2);
        
        for j=1:n
            int(j)=integral(@(y) ee(y).^n.*exp(ee(y)).*nchoosek(n,j)*(-1)^(n-j).*(ee(y)).^(-j).*gammainc(aa(y),j,'upper')*gamma(j).*y.*2./(R1^2),0,R1);
        end
        
        int_near(jj)=(1/factorial(n))*(integral(@(y) fun31(y),0,R1)+sum(int));

    end
     R_near(ii) = integral(@(y) func_ori(y),0,R1 )+1/log(2)*(sum(int_near)+integral(@(y) fun1(y),0,R1));


     %% asymptotic near 
     a_cal = a_near*beta_near*gamma_th_SIC/(a_far*beta_far-gamma_th_SIC*a_near*beta_near);
     ee=@(y) sigma_square./(a_near*beta_near.*Phy(y).*theta);
     aa=@(y) ee(y).*(a_cal+1);
     func_ori = @(y) log(a_cal+1)/log(2)*(1-(sigma_square.*gamma_th_SIC./(a_near.*beta_near.*Phy(y).*theta.*(a_index-gamma_th_SIC))).^k./k).*2.*y./R1^2;
     c=0.577215;
     fun1a= @(y) (1+ee(y)).*(-log(ee(y).*a_cal*a_index/gamma_th_SIC)-c).*y.*2./(R1^2);
       
    for jj=1:k-1
        n=jj;
        int=zeros(1,n);
        f2=@(x,y) (ee(y).*x).^n./(1+x).*exp(-ee(y).*x).*y.*2./(R1^2);
        fun31a= @(y) ee(y).^n.*(1+ee(y)).*(-1)^n.*(-log(aa(y))-c).*y.*2./(R1^2);
        
        for j=1:n
            int(j)=integral(@(y) ee(y).^n.*(1+ee(y)).*nchoosek(n,j)*(-1)^(n-j).*(ee(y)).^(-j).*gammainc(aa(y),j,'upper')*gamma(j).*y.*2./(R1^2),0,R1);
        end
        
        int_nearasy(jj)=(1/factorial(n))*(integral(@(y) fun31a(y),0,R1)+sum(int));

    end
     R_nearAsy(ii) = integral(@(y) func_ori(y),0,R1 )+1/log(2)*(sum(int_nearasy)+integral(@(y) fun1a(y),0,R1));

    %% asymptotic far
    t=@(x) sigma_square*x./(a_near*beta_near*Pt(ii)*dBR^(-alpha)*theta*a_index);
    
        
        funa=@(x) 1./(x+1)-2*t(x).^k./((R2^2-R1^2)*(x+1)*factorial(k)).*((R2^2+H^2)^(alpha*k/2+1)-(R1^2+H^2)^(alpha*k/2+1))/(alpha*k+2);
     
       
     R_farAsy(ii) = 1/log(2)*integral(@(x) funa(x),0,a_index);

end


%% plot
figure
plot(rho,R_far,'r*-' );
hold on;
plot(rho,R_farAsy,'k--' );
hold on;
ylim([0 3])
title('Far user')

figure
plot(rho,R_near,'b*-' );
hold on;
plot(rho,R_nearAsy,'k--');
hold on;
ylim([0 9])
title('Near user')
