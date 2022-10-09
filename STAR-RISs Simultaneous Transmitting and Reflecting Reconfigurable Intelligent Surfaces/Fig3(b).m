%10.1109/LCOMM.2021.3082214, Fig. 3(b)
L = 10; %L*L is the number of elements of the STAR-RIS 
N = 300;
M = 100; %M and N define the plot grid for the observation plan Sigma

table = zeros(2*N+1,2*M+1);


%x = [1,0.5,-0.5,-1,1,0.5,-0.5,-1,1,0.5,-0.5,-1,1,0.5,-0.5,-1];
%y = [1,1,1,1,0.5,0.5,0.5,0.5,-0.5,-0.5,-0.5,-0.5,-1,-1,-1,-1];
x = zeros(1,L*L);
y = zeros(1,L*L);
phi = zeros(1,L*L);
for l = 1:L*L
    x_rem = rem(l-1,L)+1;
    x(l)=x_rem/2-(L-1)/4-0.5;
    y(l)=((l-x_rem)/L+1)/2-(L-1)/4-0.5;
    phi(l) = 2*pi*y(l)/15;
    phi2(l) = -2*pi*y(l)/3;
end
%the loop evenly place each element on the STAR-RIS

%calculated the received power over each grid element (ii,jj)
for ii = 1:2*N+1
    I = (ii-N-1)/10;
    for jj = 1:2*M+1
        J = (M+1-jj)/10;
        temp = 0;
        for m = 1:L*L
            d = sqrt(x(m)^2 + (J-y(m))^2 + I^2);
            csth = abs(I/d);
            %if y(m)>0
            if rem(m,2)==0
                if I>0
                    temp = temp + (1+csth)/2*exp(1j*2*pi*d)/d*exp(1j*phi(m));
                end
            else
                if I<0
                    temp = temp + (1+csth)/2*exp(1j*2*pi*d)/d*exp(1j*phi2(m));
                end
            end           
        end
        table(ii,jj)= norm(temp);
    end
end

t2 = zeros(2*M+1,2*N+1);
for ii = 1:2*N+1
    for jj = 1:2*M+1
        t2(jj,ii)=table(ii,jj);
    end
end

%plot the radiation pattern
Z = t2;
%b = surf(Z); %use this to see 3D effect
imagesc(Z);
newmap = jet;                    %starting map
ncol = size(newmap,1);           %how big is it?
zpos = 1 + floor(2/3 * ncol);    %2/3 of way through
newmap(zpos,:) = [1 1 1];        %set that position to white
colormap(newmap);                %activate it
colorbar
%axis([0 200 0 600 0 7.3]);
xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');
%zlabel('Power Density(W/m$^2$)','Interpreter','latex');
