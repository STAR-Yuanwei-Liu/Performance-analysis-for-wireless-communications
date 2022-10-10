function [outputArg1,outputArg2]= myray2(M,D_0, flag,gamma_0)%where flag is the power ratio (eta)
it = 1e7;

gamma_1 =1e5;


gamma_db = 25:5:125; %transmit SNR
t = 10.^(gamma_db./10); %transmit SNR in dB
H = zeros(it,1);
L = size(t,2);
Pout = zeros(L,1);
Bound = zeros(L,1);
D = D_0; %D denote the quantization level where D = -1 means continuous phase-shift
thresh = 1;

coe = gamma_0^M/factorial(2*M);


for i = 1:it
    h = raylrnd(sqrt(2),M,1);
    direct = raylrnd(sqrt(2),1,1);
    if D == -1
        H(i)=sum(h);
        H(i) = H(i)/(flag^2+M)+direct*flag^2/(flag^2+M);
        continue;
    end

    theta = -pi/D + 2*pi*rand(M,1)/D;
    H(i) = sqrt((sum(h.*cos(theta)))^2+(sum(h.*sin(theta)))^2);

    
end
for j = 1:L
    Pout(j) = sum(H<sqrt(thresh*gamma_1/t(j)))/it;

    if D == -1
        if M==2
            Bound(j)= coe* (t(j)^(-M-1));
        else
            Bound(j)= coe* (t(j)^(-M-0.75));
        end
        continue;
    end
    Bound(j) = 1-exp(-thresh*gamma_0/t(j)/4/M);

end
outputArg1 = Pout;
outputArg2 = Bound;
end
