% The codes is for paper of "STAR-IOS Aided NOMA Networks: Channel Model Approximation and Performance Analysis"
%  in IEEE Transactions on Wireless Communications, vol. 21, no. 9, pp. 6861-6876, Sept. 2022, doi: 10.1109/TWC.2022.3152703.


clear 
clc
N =20; % num of RIS elements
% threhold = 70:3:100; %% no squre
threhold = 100:100:1000;
num = 5000;
Pt_dBm = 0; % transmit power dBm
P_t = 10.^((Pt_dBm-30)./10); % transmit power 

k1 = 2;
s1 = sqrt(k1./(1+k1)); % matlab paramiter one
sigma1 = sqrt( 1./ (2.*(1+k1))  ); % matlab paramiter Two

k2 = 2;
s2 = sqrt(k2./(1+k2)); % matlab paramiter one
sigma2 = sqrt( 1./ (2.*(1+k2))  ); % matlab paramiter Two






 %% simulation Rician * Rician
% for ii = 1:length(threhold)
% rician1 = random('rician',s1,sigma1,[1,num]);
% rician2 = random('rician',s2,sigma2,[1,num]);
% x = rician1.*rician2;
% CDF(ii) = sum(x<threhold(ii))/num;
% end
% plot(threhold,CDF,'r-+')
% hold on

%% simulation integrated signals
A = 1; % energy loss 0.3 rfr ; 0.7 rfl
for j = 1:length(threhold)
    tic
    th = threhold(j);
for i = 1:num
    rician1 = random('rician',s1,sigma1,[1,N]);
    rician2 = random('rician',s2,sigma2,[1,N]);
    Flag = A.*sum(rician1.*rician2).^2;
    S(i) = Flag < th;
end
 CDF(j) = sum(S)/num;
 toc
end

plot(threhold,CDF,'r*')
hold on
%% fitting
% syms t
% f=fittype('1-marcumq(sqrt(2*abs(k)*abs(u)),sqrt(2*(1+abs(k))*abs(u)).*t,*abs(u))','independent','t','coefficients',{'k','u'});  %fittype uses your self-defined function
% cfun=fit(threhold',CDF',f,'StartPoint',[N,N]); %find the value of parameters

f=fittype('(gammainc(t./b,floor(abs(a))))','independent','t','coefficients',{'a','b'});  %fittype uses your self-defined function
cfun=fit(threhold',CDF',f,'StartPoint',[N,N]); %find the value of parameters

% f=fittype('(1-marcumq(abs(v)./abs(s),t./abs(s).*t,1))','independent','t','coefficients',{'v','s'});  %fittype uses your self-defined function
% cfun=fit(threhold',CDF',f,'StartPoint',[N,N]); %find the value of parameters



out2 = cfun(threhold);
% k = 0.5;
% u = 3;
% CDF2 = 1-marcumq(sqrt(2.*k.*u),sqrt(2.*(1+k).*2).*threhold,2);
plot(threhold,out2,'b-o')
hold on
cfun