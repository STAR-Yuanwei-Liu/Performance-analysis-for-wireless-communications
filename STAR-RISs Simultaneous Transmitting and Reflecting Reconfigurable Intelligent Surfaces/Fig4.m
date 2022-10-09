%10.1109/LCOMM.2021.3082214, Fig. 4
L = 16;  %L*L is the number of elements of the STAR-RIS 
N = 2000;
d = zeros(1,2*N+1);
g = zeros(1,2*N+1);
gg = zeros(1,2*N+1);
g_0 =zeros(1,2*N+1);
%x = [1,0.5,-0.5,-1,1,0.5,-0.5,-1,1,0.5,-0.5,-1,1,0.5,-0.5,-1];
%y = [1,1,1,1,0.5,0.5,0.5,0.5,-0.5,-0.5,-0.5,-0.5,-1,-1,-1,-1];
x = zeros(1,L*L);
y = zeros(1,L*L);
phi = zeros(1,L*L);
for l = 1:L*L %the loop evenly place each element on the STAR-RIS
    x_rem = rem(l-1,L)+1;
    x(l)=x_rem/2-(L-1)/4-0.5;
    y(l)=((l-x_rem)/L+1)/2-(L-1)/4-0.5;
    phi(l) = 2*pi*y(l)/15;
end

cs = 0.5;
sn = sqrt(1-cs^2);

%calculated the received power over each position n
for n = -N:N
    d(n+N+1) = n/100;
    temp = 0;
    temp2 = 0;
    for m = 1:L*L
        dm = sqrt(x(m)^2 + (y(m)-d(n+N+1)*sn)^2 + (d(n+N+1)*cs)^2);
        csth = abs(d(n+N+1)*cs/dm);
        temp = temp+exp(1j*2*pi*dm)/dm*exp(1j*phi(m));
        temp2 =temp2 + (1+csth)/2*exp(1j*2*pi*dm)/dm*exp(1j*phi(m));
    end
    if n>0
        g(n+N+1) = norm(temp);
        gg(n+N+1) = norm(temp2);
        g_0(n+N+1) = norm(10/d(n+N+1));
    else
        g(n+N+1) = 2/3*norm(temp);
        gg(n+N+1) = 2/3*norm(temp2);
        g_0(n+N+1) = 2/3*norm(10/d(n+N+1));
    end        
end

% xx = zeros(N);
% for it = -N:N
%     xx(it+N+1) = it/100;
% end

% plot(xx(1),g(1),'b');
% hold on;
% plot(xx(1),gg(1),'r','LineWidth',2);
% hold on;
% plot(xx(1),g_0(1),'--k');

plot(d,g,'-.b');
hold on;
plot(d,gg,'r');
hold on;
plot(d,g_0,'--k');
axis([-20 20 0 6]);
xlabel('$d/\lambda$','Interpreter','latex','FontSize',16);
ylabel('Relative Channel Gain (dB)','Interpreter','latex','FontSize',16);
legend({sprintf('Near-field model \nwithout leaning factors'), sprintf('Near-field model \nwith leaning factors'),sprintf('Far-field model')})  % good
%legend({'near-field model without leaning factors','near-field model with leaning factors','far-field model'},'FontSize',12);