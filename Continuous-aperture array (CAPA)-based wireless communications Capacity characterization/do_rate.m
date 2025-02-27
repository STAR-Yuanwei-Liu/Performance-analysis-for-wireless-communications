clear;
lambda = 0.125;
SNR=10^(50/10);
d1 = 10; theta1 = pi/6; phi1 = pi/3; 
PHI1 = cos(phi1)*sin(theta1); PSI1 = sin(phi1)*sin(theta1); THETA1 = cos(theta1);
d2 = 20; theta2 = pi/6; phi2 = pi/3; 
PHI2 = cos(phi2)*sin(theta2); PSI2 = sin(phi2)*sin(theta2); THETA2 = cos(theta2);
Size = [0:0.25:5];
precision = 100;
Data = ones(6,length(Size));
for M_index = [1:1:length(Size)]
    Lx = sqrt(Size(M_index));
    Lz = sqrt(Size(M_index));
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
    rho=0;
    for j1 = 1:1:precision
        psi1=cos((2*j1-1)*pi/(2*precision));
        for j2=1:precision
            psi2=cos((2*j2-1)*pi/(2*precision));
            distance1 = sqrt((d1*PHI1 - Lx/2*psi1)^2 + (d1*PSI1 - 0)^2 + (d1*THETA1 - Lz*psi2/2)^2);
            part1 = (sqrt(d1*PSI1)/sqrt(4*pi)/(distance1^(3/2)));
            distance2 = sqrt((d2*PHI2 - Lx/2*psi1)^2 + (d2*PSI2 - 0)^2 + (d2*THETA2 - Lz*psi2/2)^2);
            part2 = (sqrt(d2*PSI2)/sqrt(4*pi)/(distance2^(3/2)));
            r_corr = r_corr + sqrt((1-psi1^2)*(1-psi2^2))*(part1*part2*cos(2*pi/lambda*(distance1 - distance2)));
            i_corr = i_corr + sqrt((1-psi1^2)*(1-psi2^2))*(part1*part2*sin(2*pi/lambda*(distance1 - distance2)));
        end
    end
    xi=0;
    p_zf=[0,0];
    if Size(M_index)>0
        rho = (pi^2*Lx*Lz/(4*precision^2))^2*(r_corr^2 + i_corr^2)/a1/a2;
        xi=(a1-a2)/(a1*a2*(1-rho)*4*pi*(120*2*pi/lambda)^2);
        p_zf=waterfill(SNR,[1/(a1*(1-rho)),1/(a2*(1-rho))]);
    end
    
    if xi>=SNR
        SNR1=SNR;
        SNR2=0;
    elseif xi<=-SNR
        SNR1=0;
        SNR2=SNR;
    else
        SNR1=(SNR-xi)/2;
        SNR2=(SNR+xi)/2;
    end

        
    C = max(log2(1 + SNR1*a1 + SNR2*a2 + SNR1*a1*SNR2*a2*(1-rho)),log2(1+SNR*a1));
    Data(1,M_index) = C;
    R1 = log2(1 + SNR1*a1);
    Data(2,M_index) = R1;
    R2 = C - R1;
    Data(3,M_index) = R2;
    R2 = log2(1 + SNR2*a2);
    Data(4,M_index) = R2;
    R1 = C - R2;
    Data(5,M_index) = R1;
    
    C_zf=log2(1+p_zf(1)*a1*(1-rho))+log2(1+p_zf(2)*a2*(1-rho));
    %C_mrc=log2(1+SNR/2*a1/(SNR/2*a2*rho+1))+log2(1+SNR/2*a2/(SNR/2*a1*rho+1));
    Data(6,M_index) = C_zf;
    %Data(7,M_index) = C_mrc;
    
end
figure
plot(Size,Data(1,:),'-*r');hold on;
plot(Size,Data(2,:),'-sb');hold on;
plot(Size,Data(3,:),'-ob');hold on;
plot(Size,Data(4,:),'-sm');hold on;
plot(Size,Data(5,:),'-om');hold on;
plot(Size,Data(6,:),'-k');hold on;