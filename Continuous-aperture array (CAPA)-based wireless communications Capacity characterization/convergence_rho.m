clear;
lambda = 0.125;
d1 = 10; theta1 = pi/6; phi1 = pi/3; 
PHI1 = cos(phi1)*sin(theta1); PSI1 = sin(phi1)*sin(theta1); THETA1 = cos(theta1);
d2 = 20; theta2 = pi/6; phi2 = pi/3; 
PHI2 = cos(phi2)*sin(theta2); PSI2 = sin(phi2)*sin(theta2); THETA2 = cos(theta2);
Size = 15;
Lx = sqrt(Size);
Lz = sqrt(Size);
precision = 1:1:25;
rho=zeros(1,length(precision));
for index = [1:1:length(precision)]
    
    a1 = 0;
    for x = [Lx/2/d1+PHI1,Lx/2/d1-PHI1]
        for z = [Lz/2/d1+THETA1,Lz/2/d1-THETA1]
            a1 = a1 + atan(x*z/PSI1/sqrt(PSI1^2+x^2+z^2))/4/pi;
        end
    end
    a2 = 0;
    for x = [Lx/2/d2+PHI2,Lx/2/d2-PHI2]
        for z = [Lz/2/d2+THETA2,Lz/2/d2-THETA2]
            a2 = a2 + atan(x*z/PSI2/sqrt(PSI2^2+x^2+z^2))/4/pi;
        end
    end
    r_corr = 0;
    i_corr = 0;
    
    for j1 = 1:1:precision(index)
        psi1=cos((2*j1-1)*pi/(2*precision(index)));
        for j2=1:precision(index)
            psi2=cos((2*j2-1)*pi/(2*precision(index)));
            distance1 = sqrt((d1*PHI1 - Lx/2*psi1)^2 + (d1*PSI1 - 0)^2 + (d1*THETA1 - Lz*psi2/2)^2);
            part1 = (sqrt(d1*PSI1)/sqrt(4*pi)/(distance1^(3/2)));
            distance2 = sqrt((d2*PHI2 - Lx/2*psi1)^2 + (d2*PSI2 - 0)^2 + (d2*THETA2 - Lz*psi2/2)^2);
            part2 = (sqrt(d2*PSI2)/sqrt(4*pi)/(distance2^(3/2)));
            r_corr = r_corr + sqrt((1-psi1^2)*(1-psi2^2))*(part1*part2*cos(2*pi/lambda*(distance1 - distance2)));
            i_corr = i_corr + sqrt((1-psi1^2)*(1-psi2^2))*(part1*part2*sin(2*pi/lambda*(distance1 - distance2)));
        end
    end
    
    rho(index) =((pi^2*Size/(4*precision(index)^2))^2*(r_corr^2 + i_corr^2)/a1/a2);
 
  
end
figure
plot(precision,rho,'-*r');hold on;
xlim([0,25]);
ylim([0,1]);