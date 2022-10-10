[a,bb]=myray2(3,-1,0.1,5e5); %diversity-r
semilogy(a,"ok");
hold on;
[a,cc]=myray2(3,-1,0.5,5e5); %upper (primary)
semilogy(a,"sk");
[a,dd]=myray2(3,-1,1,7.5e5); %upper (primary)
semilogy(a,"dk");
hold on;



%%%%%%%%%%%%%%
[a,bb2]=myray2(2,-1,0.1,5e5); %diversity-r
semilogy(a,"ob");

hold on;
[a,cc]=myray2(2,-1,0.5,5e5); %upper (primary)
semilogy(a,"sb");


[a,dd2]=myray2(2,-1,1,5e5); %upper (primary)
semilogy(a,"db");
hold on;

semilogy(bb,"k");

semilogy(bb2,"b");



axis([3 8 0.000001 1]);
xlabel('SNR (dB)','Interpreter','latex');
ylabel('Outage Probability','Interpreter','latex');