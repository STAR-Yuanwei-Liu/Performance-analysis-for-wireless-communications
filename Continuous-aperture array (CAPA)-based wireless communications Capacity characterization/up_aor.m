clear;
lambda = 0.125;
SNR1 = 10^(30/10);
SNR2 = 10^(40/10);
d1 = 10; theta1 = pi/6; phi1 = pi/3; 
PHI1 = cos(phi1)*sin(theta1); PSI1 = sin(phi1)*sin(theta1); THETA1 = cos(theta1);
d2 = 20; theta2 = pi/6; phi2 = pi/3; 
PHI2 = cos(phi2)*sin(theta2); PSI2 = sin(phi2)*sin(theta2); THETA2 = cos(theta2);
precision = 30;
AOR_Term = [0:0.05:1];
Data = zeros(6,length(AOR_Term));
for AOR_Index = [1:1:length(AOR_Term)]
    AOR = AOR_Term(AOR_Index);
    Size = 1;
    Lx = sqrt(Size);
    Lz = sqrt(Size);
    %% SPD 
    a1 = 0;
    for x = [Lx/2/d1+PHI1,Lx/2/d1-PHI1]
        for z = [Lz/2/d1+THETA1,Lz/2/d1-THETA1]
            a1 = a1 + AOR*atan(x*z/PSI1/sqrt(PSI1^2+x^2+z^2))/4/pi;
        end
    end
    a2 = 0;
    for x = [Lx/2/d2+PHI2,Lx/2/d2-PHI2]
        for z = [Lz/2/d2+THETA2,Lz/2/d2-THETA2]
            a2 = a2 + AOR*atan(x*z/PSI2/sqrt(PSI2^2+x^2+z^2))/4/pi;
        end
    end

    [x, w] = GaussLegendre(precision);
    r_corr = 0;
    i_corr = 0;
    for index = [1:1:precision]
        distance1 = sqrt((d1/d2*PHI1 - Lx/d2/2*x(index))^2 + (d1/d2*PSI1 - 0)^2 + (d1/d2*THETA1 - Lz/d2*x/2).^2);
        part1 = (sqrt(d1/d2*PSI1)/sqrt(4*pi)./(distance1.^(3/2)));
        distance2 = sqrt((PHI2 - Lx/d2/2*x(index))^2 + (PSI2 - 0)^2 + (THETA2 - Lz/d2*x/2).^2);
        part2 = (sqrt(PSI2)/sqrt(4*pi)./(distance2.^(3/2)));
        r_corr = r_corr + Lx/2/d2*w(index)*(Lz/d2/2*transpose(w)*(part1.*part2.*cos(2*pi/lambda*d2*(distance1 - distance2))));
        i_corr = i_corr + Lx/2/d2*w(index)*(Lz/d2/2*transpose(w)*(part1.*part2.*sin(2*pi/lambda*d2*(distance1 - distance2))));
    end
    rho = AOR^2*(r_corr^2 + i_corr^2)/a1/a2;
    C = log2(1 + SNR1*a1 + SNR2*a2 + SNR1*a1*SNR2*a2*(1-rho));
    if AOR>0 
        Data(1,AOR_Index) = C;
    end
    %% CAPA (Upper Limit)
    Lx = sqrt(Size);
    Lz = sqrt(Size);
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
    rho = (pi^2*Lx*Lz/(4*precision^2))^2*(r_corr^2 + i_corr^2)/a1/a2;
    C = log2(1 + SNR1*a1 + SNR2*a2 + SNR1*a1*SNR2*a2*(1-rho));
    Data(2,AOR_Index) = C;
    
    Size = 3;
    Lx = sqrt(Size);
    Lz = sqrt(Size);
    %% SPD 
    a1 = 0;
    for x = [Lx/2/d1+PHI1,Lx/2/d1-PHI1]
        for z = [Lz/2/d1+THETA1,Lz/2/d1-THETA1]
            a1 = a1 + AOR*atan(x*z/PSI1/sqrt(PSI1^2+x^2+z^2))/4/pi;
        end
    end
    a2 = 0;
    for x = [Lx/2/d2+PHI2,Lx/2/d2-PHI2]
        for z = [Lz/2/d2+THETA2,Lz/2/d2-THETA2]
            a2 = a2 + AOR*atan(x*z/PSI2/sqrt(PSI2^2+x^2+z^2))/4/pi;
        end
    end

    [x, w] = GaussLegendre(precision);
    r_corr = 0;
    i_corr = 0;
    for index = [1:1:precision]
        distance1 = sqrt((d1/d2*PHI1 - Lx/d2/2*x(index))^2 + (d1/d2*PSI1 - 0)^2 + (d1/d2*THETA1 - Lz/d2*x/2).^2);
        part1 = (sqrt(d1/d2*PSI1)/sqrt(4*pi)./(distance1.^(3/2)));
        distance2 = sqrt((PHI2 - Lx/d2/2*x(index))^2 + (PSI2 - 0)^2 + (THETA2 - Lz/d2*x/2).^2);
        part2 = (sqrt(PSI2)/sqrt(4*pi)./(distance2.^(3/2)));
        r_corr = r_corr + Lx/2/d2*w(index)*(Lz/d2/2*transpose(w)*(part1.*part2.*cos(2*pi/lambda*d2*(distance1 - distance2))));
        i_corr = i_corr + Lx/2/d2*w(index)*(Lz/d2/2*transpose(w)*(part1.*part2.*sin(2*pi/lambda*d2*(distance1 - distance2))));
    end
    rho = AOR^2*(r_corr^2 + i_corr^2)/a1/a2;
    C = log2(1 + SNR1*a1 + SNR2*a2 + SNR1*a1*SNR2*a2*(1-rho));
    %Data(3,AOR_Index) = C;
    if AOR>0 
        Data(3,AOR_Index) = C;
    end
    %% CAPA (Upper Limit)
    Lx = sqrt(Size);
    Lz = sqrt(Size);
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
    rho = (pi^2*Lx*Lz/(4*precision^2))^2*(r_corr^2 + i_corr^2)/a1/a2;
    C = log2(1 + SNR1*a1 + SNR2*a2 + SNR1*a1*SNR2*a2*(1-rho));
    Data(4,AOR_Index) = C;
    
    Size = 5;
    Lx = sqrt(Size);
    Lz = sqrt(Size);
    %% SPD 
    a1 = 0;
    for x = [Lx/2/d1+PHI1,Lx/2/d1-PHI1]
        for z = [Lz/2/d1+THETA1,Lz/2/d1-THETA1]
            a1 = a1 + AOR*atan(x*z/PSI1/sqrt(PSI1^2+x^2+z^2))/4/pi;
        end
    end
    a2 = 0;
    for x = [Lx/2/d2+PHI2,Lx/2/d2-PHI2]
        for z = [Lz/2/d2+THETA2,Lz/2/d2-THETA2]
            a2 = a2 + AOR*atan(x*z/PSI2/sqrt(PSI2^2+x^2+z^2))/4/pi;
        end
    end

    [x, w] = GaussLegendre(precision);
    r_corr = 0;
    i_corr = 0;
    for index = [1:1:precision]
        distance1 = sqrt((d1/d2*PHI1 - Lx/d2/2*x(index))^2 + (d1/d2*PSI1 - 0)^2 + (d1/d2*THETA1 - Lz/d2*x/2).^2);
        part1 = (sqrt(d1/d2*PSI1)/sqrt(4*pi)./(distance1.^(3/2)));
        distance2 = sqrt((PHI2 - Lx/d2/2*x(index))^2 + (PSI2 - 0)^2 + (THETA2 - Lz/d2*x/2).^2);
        part2 = (sqrt(PSI2)/sqrt(4*pi)./(distance2.^(3/2)));
        r_corr = r_corr + Lx/2/d2*w(index)*(Lz/d2/2*transpose(w)*(part1.*part2.*cos(2*pi/lambda*d2*(distance1 - distance2))));
        i_corr = i_corr + Lx/2/d2*w(index)*(Lz/d2/2*transpose(w)*(part1.*part2.*sin(2*pi/lambda*d2*(distance1 - distance2))));
    end
    rho = AOR^2*(r_corr^2 + i_corr^2)/a1/a2;
    C = log2(1 + SNR1*a1 + SNR2*a2 + SNR1*a1*SNR2*a2*(1-rho));
    if AOR>0 
        Data(5,AOR_Index) = C;
    end
    %% CAPA (Upper Limit)
    Lx = sqrt(Size);
    Lz = sqrt(Size);
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
    rho = (pi^2*Lx*Lz/(4*precision^2))^2*(r_corr^2 + i_corr^2)/a1/a2;
    C = log2(1 + SNR1*a1 + SNR2*a2 + SNR1*a1*SNR2*a2*(1-rho));
    Data(6,AOR_Index) = C;
end
figure
plot(AOR_Term,Data(1,:),'-sr');hold on;
plot(AOR_Term,Data(2,:),'-ob');hold on;
plot(AOR_Term,Data(3,:),'--sb');hold on;
plot(AOR_Term,Data(4,:),'--om');hold on;
plot(AOR_Term,Data(5,:),':sm');hold on;
plot(AOR_Term,Data(6,:),':om');hold on;