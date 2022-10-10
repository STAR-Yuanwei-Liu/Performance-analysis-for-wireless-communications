function [outputArg1,outputArg2]= myray(M,D_0, flag,gamma_0)
it = 1e5;
gamma_db = 25:3.5:105;
t = 10.^(gamma_db./10);
H = zeros(it,1);
L = size(t,2);
Pout = zeros(L,1);
Bound = zeros(L,1);
D = D_0;
thresh = 1;

Rician_factor = 1;
pd_1 = makedist('Rician','s',sqrt(Rician_factor/(Rician_factor+1)),'sigma',sqrt(0.5/(Rician_factor+1)));
pd_2 = makedist('Rician','s',sqrt(2/3),'sigma',sqrt(0.5/3));



coe = gamma_0^M/factorial(2*M)/6;

for i = 1:it
    h = raylrnd(sqrt(2),M,1);

    if D == -1
        H(i)=sum(h);
        continue;
    end

    theta = -pi/D + 2*pi*rand(M,1)/D;
    H(i) = sqrt((sum(h.*cos(theta)))^2+(sum(h.*sin(theta)))^2);

end
for j = 1:L
    Pout(j) = sum(H<sqrt(thresh*gamma_0/t(j)))/it;
    if D == -1
        Bound(j)= coe* (t(j)^(-M));
        continue;
    end
    Bound(j) = 1-exp(-thresh*gamma_0/t(j)/4/M);
end
outputArg1 = Pout;
outputArg2 = Bound;
end
