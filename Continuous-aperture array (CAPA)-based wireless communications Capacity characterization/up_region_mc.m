clear;

Lx=1;
Lz=1;

lambda=0.125;

d=lambda/3;
A=d^2;
A=lambda^2/12;

Mx=round(Lx/d);
Mz=round(Lz/d);

SNR1 = 10^(30/10);
SNR2 = 10^(40/10);

theta_b=pi/3;
theta_e=pi/3;
phi_b=pi/6;
phi_e=pi/6;

Phi_b=sin(phi_b)*cos(theta_b);
Phi_e=sin(phi_e)*cos(theta_e);
Psi_b=sin(theta_b)*sin(phi_b);
Psi_e=sin(theta_e)*sin(phi_e);
Omega_b=cos(phi_b);
Omega_e=cos(phi_e);

r1=10;
r2=20;
p_1=[r1*sin(phi_b)*cos(theta_b),r1*sin(theta_b)*sin(phi_b),r1*cos(phi_b)];
p_2=[r2*sin(phi_e)*cos(theta_e),r2*sin(theta_e)*sin(phi_e),r2*cos(phi_e)];


hb_mat=zeros(Mx,Mz);
he_mat=zeros(Mx,Mz);

for i=1:Mx
    for j=1:Mz
        p_ij=[(i-(Mx+1)/2)*d,0,(j-(Mz+1)/2)*d];
        n_b=norm(p_1-p_ij);
        n_e=norm(p_2-p_ij);
        hb_mat(i,j)=sqrt(A*r1*sin(theta_b)*sin(phi_b)/(4*pi*n_b^3))*exp(-1i*2*pi/lambda*n_b);
        he_mat(i,j)=sqrt(A*r2*sin(theta_e)*sin(phi_e)/(4*pi*n_e^3))*exp(-1i*2*pi/lambda*n_e);

    end
end
h1=reshape(hb_mat,[],1);
h2=reshape(he_mat,[],1);
g1=norm(h1)^2;
g2=norm(h2)^2;
rho=abs(h1'*h2)^2/(g1*g2);

C = log2(1 + SNR1*g1 + SNR2*g2 + SNR1*g1*SNR2*g2*(1-rho));

R21_1 = log2(1 + SNR1*g1);

R21_2= C - R21_1;

R12_2 = log2(1 + SNR2*g2);

R12_1 = C - R12_2;

figure
plot([0,R12_1,R21_1,R21_1],[R12_2,R12_2,R21_2,0]);
hold on;

Z=ones(Mx*Mz,Mx*Mz);
for i=1:Mx*Mz
    for j=1:Mx*Mz
        if i~=j
            [i,j]
            x1=mod(i,Mx);
            y1=ceil(i/Mz);
            x2=mod(j,Mx);
            y2=ceil(j/Mz);
            delta=sqrt((x1-x2)^2+(y1-y2)^2)*d;
            Z(i,j)=0.1*exp(-1j*2*pi/lambda*delta)/delta^2;
        end

    end
end
C_matrix=100*(Z+50*eye(Mx*Mz))^(-1);
h1=C_matrix'*h1;
h2=C_matrix'*h2;
g1=norm(h1)^2;
g2=norm(h2)^2;
rho=abs(h1'*h2)^2/(g1*g2);
C = log2(1 + SNR1*g1 + SNR2*g2 + SNR1*g1*SNR2*g2*(1-rho));

R21_1 = log2(1 + SNR1*g1);

R21_2= C - R21_1;

R12_2 = log2(1 + SNR2*g2);

R12_1 = C - R12_2;

plot([0,R12_1,R21_1,R21_1],[R12_2,R12_2,R21_2,0]);
hold on;