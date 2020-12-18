function [ n, W,SumP, quotientA,w] = main( sigmaA,sigmaF,sigmaS,nNeutron,tMax)
    % nem analóg szimuláció szórásminimalizálás
warning('off','all')
close all

nu = 2.5; % hasadási sokszorozódás
sigmaT = sigmaA + sigmaS + sigmaF;

pA=sigmaA/sigmaT;  %valségek [0,1]
pF=sigmaF/sigmaT;
pS=sigmaS/sigmaT;

SumP = zeros(1,200);
W= zeros(1,nNeutron);
for i = 1:200 % i a q-kra megy
    W= zeros(1,nNeutron); % hasadások száma (fi)
    w = ones(1,nNeutron); % neutronok súlyai (populációkontroll)
    nPop = zeros(1,tMax); % neutronpopuláció az idõben
    var = zeros(1,tMax); % relatív szórások vektora (pop kontroll)
    varW = zeros(1,tMax); % hasadások szórása (fi)
    %q(i) = 0 + 0.005*i; %teljes térkép
    quotientA(i) = pF*nu + pS;
     %q(i) = 0.5+0.0005*i; % pop kontroll közeli
     q(i) = 0.5+0.0005*i; % fi közeli
    for t = 1:tMax %idõléptetés
        for j = 1:nNeutron %neutronléptetés
            w(j) = (1 - pA)*w(j); %túlélés valsége
            if rand  < q(i) % hasadás, most 2.5 neutron keletkezik egy hasadásban
                w(j) = 1/q(i)*pF/(pF + pS)*w(j);
                SumP(i) = SumP(i) + w(j);
                W(j) = W(j) + w(j);
                w(j) = nu*w(j);
            else
                w(j) = w(j)*(1 - pF/(pF + pS))/(1 - q(i)); % szóródik, súlykorrekció itt is kell
            end
        end
        nPop(t) = sum(w); % neutronpopuláció a t. idõben
        var(t) = sum(w.^2)/sum(w)^2 - 1/nNeutron; %súlyok szórása
        varW(t) = sum(W.^2)/sum(W)^2 - 1/nNeutron; %súlyok szórása
        Ew(t) = sum(w)/nNeutron; %súlyok várható értéke, átlagos súly
    end
    n(i) = nPop(tMax); % neutronok száma a szimuláció végén
    qvar(i) = var(tMax); %szimuláció végén a súlyok szórása
    qvarW(i) = varW(tMax); %szimuláció végén a súlyok szórása
    quotientS(i) = nthroot(n(i)/nNeutron,tMax); % szimuláció kvóciense, szimuláció által adott átlag szaporulat
    
end
figure(1)
plot(q,qvar,'b.-') %súlyok szórásnégyzete
%plot(q,sqrt(qvar),'b.-') %súlyszórás
grid on
xlabel('q','fontsize',18)
ylabel('Relative variance (RSD^2)','fontsize',18)
title('RSD^2 - q','fontsize',18)
legend('Relative variance of the neutron weights','fontsize',18)

figure(2)
%plot(q,qvarW,'r.-') %súlyok szórásnégyzete
plot(q,sqrt(qvarW),'r.-') %súlyszórás
grid on
xlabel('q','fontsize',18)
ylabel('Relative variance (RSD^2)','fontsize',18)
title('RSD^2 - q','fontsize',18)
legend('Relative variance of the power output','fontsize',18)

figure(3)
hold on
plot(q,quotientS,'g.-') % szimuláció alapján számolt átlagkvóciens 
plot(q,quotientA,'r') % analitikus kvóciens (várható értéke egy lépés után)
grid on
legend('Simulated quotient','Analytic quotient','fontsize',18)
xlabel('q','fontsize',18)
ylabel('Average population growth in one step','fontsize',14)
title('Quotient - q','fontsize',18)
hold off 

figure(4)
hold on
plot(q,n,'m.-') % végsõ neutronszám
plot(q,nNeutron*quotientA.^tMax,'r'); %elméleti végsõ neutronszám
grid on
legend('Simulated final population','Analytic final population','fontsize',18)
xlabel('q','fontsize',18)
ylabel('Final population','fontsize',18)
title('Final popultion - q','fontsize',18)
hold off


end

