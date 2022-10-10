theta = 0:0.001:2*pi;
len = size(theta,2);
M = 12;% change M: number of element
R = 100;%*lambda
L = (M-1)/2;
y = -L/2:0.5:L/2;
rho = zeros(1,len);
sn = sind(30);

for it = 1:len
    for m = 1:M
        dm = (R*cos(theta(it)))^2+(R*sin(theta(it))-y(m))^2;% distance between user and element m
        dm = sqrt(dm);
        con = mod(pi*sn*m,2*pi);
        pha = mod(2*pi*dm,2*pi);
        phase_error =  (rand()-0.5)*pi*2; %change phase_error for the corresponding DP-PSC and PS-PSC strategies
        leaning_factor = (1+abs(cos(theta(it))))/2;
        if (theta(it)<pi/2 || theta(it)>3*pi/2)
            %uncomment for ploting the random PSC:
            %rho(it) = rho(it) + leaning_factor*exp(1j*phase_error)*exp(2j*pi*dm)/dm;
            rho(it) = rho(it) + leaning_factor*exp(2j*pi*dm)/dm;%TR-PSC
        else
            rho(it) = rho(it) + leaning_factor*exp(1j*pi*sn*m)*exp(2j*pi*dm)/dm;%TR-PSC
            %rho(it) = rho(it) + leaning_factor*exp(1j*phase_error)*exp(1j*pi*sn*m)*exp(2j*pi*dm)/dm;
        end
    end
    %rho(it) = norm(rho(it)) * (1+abs(cos(theta(it))))/2;
    rho(it) = norm(rho(it));
end


%rho = sin(2*theta).*cos(2*theta);
polarplot(theta,rho);
rlim([0 0.175])