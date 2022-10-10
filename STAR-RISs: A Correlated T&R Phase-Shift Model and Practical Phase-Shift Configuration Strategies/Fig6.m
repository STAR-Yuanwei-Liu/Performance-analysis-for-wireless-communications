%Fig6
[a,bb]=myray(4,4,1,1e5); %diversity-r
semilogy(a,"or");
hold on;
[a,cc]=myray(4,-1,1,1e5); %upper (primary)

hold on;

[a,bb]=myray(4,4,1,5e5); %diversity-t
semilogy(a,"sr");
[a,ee]=myray(4,-1,1,5e5);

[a2,bb2]=myray(4,-1,1,3e5); %de-r
semilogy(a2,">b");

[a2,cc2]=myray(4,-1,1,8e5);%de-t
semilogy(a2,"<b");


[a2,dd2]=myray(4,1,1,1e5);%random-r
semilogy(a2,"^",'MarkerEdgeColor','g','MarkerFaceColor',[0 0.7 0.7]);
[a2,ee2]=myray(4,1,1,5e5);%random-t
semilogy(a2,"v",'MarkerEdgeColor','g','MarkerFaceColor',[0 0.7 0.7]);

semilogy(bb,"k");
semilogy(dd2,"k");

semilogy(bb2,"b");
semilogy(cc,"k");
semilogy(ee,"k");



semilogy(cc2,"b");

semilogy(ee2,"k");


axis([1 9 0.0001 1]);
xlabel('SNR (dB)','Interpreter','latex');
ylabel('Outage Probability','Interpreter','latex');
