warning('off','all')
close all
clear all
pF = 0.25;
pS = 0.3;
pA = 1 - pF - pS;
betai=[21 142 128 257 75 27]*1e-5;
T12i=[55.7 22.7 6.2 2.3 0.615 0.23];
beta=sum(betai);
lambda=1/( sum(T12i.*betai)/beta );
for m = 1:6
nu = 2.5; % hasadási sokszorozódás
nNeutron = 10000;
tMax = 15*m;
resolution = 5;
iMax = 10; % hogy átlagolni tudjon iMaxszor lefut
maxStep = 20; % eddig vizsgálja hogy ennyi lépésbõl kiugrik-e
qopt = 1 - pS;
betaopt = 1 - pF*nu*(1 - beta)/qopt;
lambdaopt = lambda*pF*nu*beta/(qopt*betaopt);
quotientprompt = pS + (1 - beta)*pF*nu;
quotientall = pS + pF*nu;
if rem(tMax,resolution) == 0
    interval = resolution:resolution:tMax; % idõintervallumok
else %hibaüzenet
    error('A max idõ legyen osztható a felbontással!');
    exit(0);
end
JumpNChain = zeros(iMax,maxStep);
for i = 1:iMax
SumP = zeros(length(interval),nNeutron);
SumPossz = 0;
nStep = zeros(maxStep,1);
lefut = 0;
        wT=zeros(length(interval),nNeutron);
        w = ones(1,nNeutron); % neutronok súlyai
        W = zeros(length(interval),nNeutron); % súlyok bankja
        for j = 1:nNeutron %neutronléptetés
            step = 0; %hányat ugrik mire kirepül
            t = 0;
           k = 1;
            while t < max(interval)
              rr=rand;
              w(j) = (1 - pA)*w(j);
                if rr  < pF/(1-pA) % hasadás,prompt neutronok csak
                   SumP(k,j) = SumP(k,j) + w(j); % leadott P részecskére, intervallumra bontva
                    SumPossz = SumPossz + w(j); % leadott összteljesítmény, egész, összes
                    W(k,j) = W(k,j) + w(j); % leadott teljesítmény neutrononként, intervallumonként eltárolva: sum(W(k,:)) = SumP(k) teljesül
                    w(j) = w(j)*nu;
                    r = rand;
                    if r < beta % itt keletkezik késsõneutron: ez csak idõben ugrik
                        dt=-1/lambda*log(rand); % életidõ sorsolás torzítottan
                        step = step + 1;
                        for n = 1:maxStep % meddig vizsgálja
                        if dt > (tMax - t) && step == n
                        nStep(n) = nStep(n) + sum(SumP(:,j));
                        JumpNChain(i,n) = JumpNChain(i,n) + 1;
                        end
                        end
                        t = t + dt;
                        if t > interval(k) % ha az idõ túllépi az adott intervallumot, akkor a következõbe lép
                            t = interval(k); % szakaszhatárra állítás
                            wT(k,j) = w(j); % az adott intervallum határON a neutronlánc súlya
                            if t < max(interval)  %!!!!!! <= nem kell csak <, mert akkor az utolsónál is belép, de ott nem kell
                                dt = -1/lambda*log(rand);
                                t = t + dt;   
                            end
                            if k + 1 + floor(dt/resolution) <= length(interval)%ha nem ugrik ki
                                k = k + 1 + floor(dt/resolution); % egynél többet lép itt esetenként! nem 5tel hanem a lépéssel kell leosztani
                            else % (10,20),30,40 re is ugorhat, ilyenkor 20ra kell állítani, különben k index túllép
                                k = length(interval);
                            end
                        elseif  t < interval(k) % ha nem lépi túl az intervallumhatárt, akkor is kell korrekció 
                        end
                    else
                        lefut = lefut + 1;
                    end
                end
            end
        end
sumPower(i) = SumPossz;
for n = 1:maxStep
nStepNorm(i,n) = nStep(n)/SumPossz;
sumNJumpPower(i,n) = nStep(n);
end
end
for n = 1:maxStep
AvgNJumpPowerNormed(n) = sum(nStepNorm(:,n))/iMax; % MENNYI TELJESÍTMÉNY NORMÁLVA
AvgNJumpChainNormed(n) = sum(JumpNChain(:,n))/iMax/nNeutron; % HÁNY LÁNC NORMÁLVA
AvgNJumpPower(n) = sum(sumNJumpPower(:,n))/iMax; % MENNYI TELJESÍTMÉNY NEM NORMÁLVA
end
AvgSumP = sum(sumPower)/iMax;
analyticalQuotient = pF*nu*(1 - beta) + pS;
analyticalSumP = nNeutron*pF/(1 - (pS + nu*pF*(1 - beta)));
for n = 1:maxStep
analyticalValues(n) = (lambda*tMax)^(n - 1)/factorial(n - 1)*exp(-lambda*tMax);
end
 figure(1)
 hold on
 plot(1:maxStep,AvgNJumpPowerNormed,'MarkerSize',10,'LineWidth',1)
 %plot(1:maxStep,AvgNJumpChainNormed,'MarkerSize',10,'LineWidth',1)
 xlabel('Escapes with "n" steps','fontsize',18)
 grid on
 ylabel('Proportion of the power','fontsize',16)
legend('Interval length = 15','Interval length = 30','Interval length = 45','Interval length = 60','Interval length = 75','Interval length = 90');
 title({'Probability of escaping the system with','"n" steps represented by','normalized quantities'},'fontsize',12)
 hold off

figure(2)
subplot(2,3,m)
hold on
plot(1:maxStep,AvgNJumpPowerNormed,'r','MarkerSize',10,'LineWidth',3)
plot(1:maxStep,AvgNJumpChainNormed,'b-.','MarkerSize',10,'LineWidth',3)
plot(1:maxStep,analyticalValues,'g','MarkerSize',10,'LineWidth',1)
str = sprintf('Interval length = %d', tMax);
title(str,'FontSize',14);
xlabel('Escapes with "n" steps','fontsize',14)
ylabel('Proportion of the power/chains','fontsize',15)
grid on
end
legend('Power output distribution','Chain distribution','Analytical proportions')
set(legend,'FontSize',16);
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
 text(0.5, 1,'\bf Probability of escaping the system with "n" steps represented by normalized proportions','HorizontalAlignment','center','VerticalAlignment', 'top','fontsize',16)

hold off