function [ n, nFission, quotientA] = main( sigmaA,sigmaF,sigmaS,nNeutron,tMax)
    % nem analóg szimuláció szórásminimalizálás
warning('off','all')
close all

nu = 2.5; % hasadási sokszorozódás
sigmaT = sigmaA + sigmaS + sigmaF;

pA=sigmaA/sigmaT;  %valségek [0,1]
pF=sigmaF/sigmaT;
pS=sigmaS/sigmaT;
for i = 1:200
    w = ones(1,nNeutron); % neutronok súlyai
    f = zeros(1,tMax); % hasadások száma
    nPop = zeros(1,tMax); % neutronpopuláció az idõben
    var = zeros(1,tMax); % relatív szórások vektora
    q(i) = 0 + 0.005*i;
    quotientA(i) = pF*nu + pS;
    %q(i) = 0.5+0.0005*i;
    for t = 1:tMax %idõléptetés
        for j = 1:nNeutron %neutronléptetés
            w(j) = (1 - pA)*w(j); %túlélés valsége
            if rand  < q(i) % hasadás, most 2.5 neutron keletkezik egy hasadásban
                w(j) = 1/q(i)*pF/(pF + pS)*nu*w(j);
                f(t) = f(t) + 1;
            else
                w(j) = w(j)*(1 - pF/(pF + pS))/(1 - q(i)); % szóródik, súlykorrekció itt is kell
            end
        end
        nPop(t) = sum(w); % neutronpopuláció a t. idõben
        var(t) = sum(w.^2)/sum(w)^2 - 1/nNeutron; %súlyok szórása
        Ew(t) = sum(w)/nNeutron; %súlyok várható értéke, átlagos súly
    end
    n(i) = nPop(tMax); % neutronok száma a szimuláció végén
    nFission(i) = sum(f); % hasadások száma
    qvar(i) = var(tMax); %szimuláció végén a súlyok szórása
    %varP(i) = sum(nPop.^2)/sum(nPop)^2 - 1/nNeutron; % neutronpopuláció
    quotientS(i) = nthroot(n(i)/nNeutron,tMax); % szimuláció kvóciense, szimuláció által adott átlag szaporulat
    
end
figure(1)
plot(q,qvar,'b.-') %súlyok szórásnégyzete
%plot(q,sqrt(qvar),'b.-') %súlyszórás
grid on
xlabel('q')
ylabel('Relative deviation')
title('Relative deviation - q')
legend('Relative variance of the test')


figure(2)
hold on
plot(q,nFission,'k.-') % összes hasadásszám egy szimulációban
plot(q,q*nNeutron*tMax,'r') %összes hasadásszám analitikus értéke
grid on
xlabel('q')
ylabel('Sum of fissions')
title('Sum of fissions - q')
hold off

figure(3)
hold on
plot(q,quotientS,'g.-') % szimuláció alapján számolt átlagkvóciens 
plot(q,quotientA,'r') % analitikus kvóciens (várható értéke egy lépés után)
grid on
legend('Simulated quotient','Analytic quotient')
xlabel('q')
ylabel(' Expected population growth in one step (quotient)')
title('Quotient - q')
hold off 

figure(4)
hold on
plot(q,n,'m.-') % végsõ neutronszám
plot(q,nNeutron*quotientA.^tMax,'r'); %elméleti végsõ neutronszám
grid on
legend('Simulated final population','Analytic final population')
xlabel('q')
ylabel('Final population')
title('Final popultion - q')
hold off


end

