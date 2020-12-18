for j = 1:8
delete(figure(j));
end
warning('off','all')
preSample = [1000,10000,25000,50000,75000,100000,150000,200000,300000,400000,500000,750000,1000000,1500000,2000000];
%preSample = [1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,11000,12000,13000,14000,15000];
%for k = 1:500
%preSample(k) = k*25000;
%end

for task = 1:4
    check = 0;
%task = 4;
for i = 1:15
    nSample(i) = preSample(i);
[integral(task,i),monteCarloIntegral(task,i),relDev(task,i),time(task,i),check] = MCintegral(0.6,1.5,nSample(i),1,task,i,check);
fom(task,i) = 1/(relDev(task,i)^2*time(task,i)); % figure of merit
end
end
figure(5)
semilogx(nSample,relDev,'o-','MarkerSize',10,'LineWidth',2)
legend('Task 1','Task 2','Task 3','Task 4')
xlabel('Number of samples','fontsize',18)
ylabel('Relative deviation','fontsize',18)
title('Relative deviation~Number of samples','fontsize',18)
set(legend,'FontSize',14);
grid on

figure(6)
semilogx(nSample,monteCarloIntegral,'o-',nSample,integral,'g-.','MarkerSize',10,'LineWidth',2)
legend('MC integral task 1','MC integral task 2','MC integral task 3','MC integral task 4','Analytic integral')
xlabel('Number of samples','fontsize',18)
ylabel('Integral result','fontsize',18)
title('Integral approximation~Number of samples','fontsize',18)
set(legend,'FontSize',14);
grid on

figure(7)
plot(nSample,fom,'o-','MarkerSize',10,'LineWidth',2)
legend('FOM task 1','FOM task 2','FOM task 3','FOM task 4')
xlabel('Number of samples','fontsize',18)
ylabel('Figure of merit','fontsize',18)
title('FOM~Number of samples','fontsize',18)
set(legend,'FontSize',14);
grid on

figure(8)
plot(nSample,time,'o-','MarkerSize',10,'LineWidth',2)
legend('Run time task 1','Run time task 2','Run time task 3','Run time task 4')
xlabel('Number of samples','fontsize',18)
ylabel('Figure of merit','fontsize',18)
title('Run time~Number of samples','fontsize',18)
set(legend,'FontSize',14);
grid on

figure(9)
semilogx(nSample,fom,'o-','MarkerSize',10,'LineWidth',2)
legend('FOM task 1','FOM task 2','FOM task 3','FOM task 4')
xlabel('Number of samples','fontsize',18)
ylabel('Figure of merit','fontsize',18)
title('FOM~Number of samples','fontsize',18)
set(legend,'FontSize',14);
grid on

figure(10)
semilogx(nSample,time,'o-','MarkerSize',10,'LineWidth',2)
legend('Run time task 1','Run time task 2','Run time task 3','Run time task 4')
xlabel('Number of samples','fontsize',18)
ylabel('Figure of merit','fontsize',18)
title('Run time~Number of samples','fontsize',18)
set(legend,'FontSize',14);
grid on