function [ n, nFission,quotientA,quotientS,deviationA,deviationS] = neutronFreePath( sigmaA,sigmaF,sigmaS,nNeutron,tMax,fP,method )
%tMax - maximum idõléptetés
%fP = 0 - nincs szabadúthossz sorsolás 
%fP = 1 - van
% method = 0 - mindhárom módszer lefut
% method = 1/2/3 - 1./2./3. módszer fut csak le

warning('off','all')
close all
if (fP == 0 || fP == 1) && (method == 0 || method == 1 || method == 2 || method == 3) && nNeutron > 0 && tMax > 0 && sigmaA >= 0 && sigmaF >= 0 && sigmaS >= 0 
for m = 1:3
    nu = 2.5; % hasadási sokszorozódás
    sigmaT = sigmaA + sigmaS + sigmaF; % teljes hatáskeresztmetszet
    w = ones(1,nNeutron); % neutronok súlyai (1,1,1,....) kezdetben
    f = zeros(1,tMax); % hasadások száma adott lépésben
    nPop = zeros(1,tMax); % neutronpopuláció/teljesítmény az idõben
    var = zeros(1,tMax); % súlyok relatív szórásának vektora
    pA=sigmaA/sigmaT;  %valségek [0,1]
    pF=sigmaF/sigmaT;
    pS=sigmaS/sigmaT;
    freePath = zeros(1,tMax); % szabadúthossz (idõ most????)
    %q = (nu*pF)/(pS + nu*pF); %ezzel a súlyszórás minimalizálódik, népességszabályozás
    q = (sigmaT-sigmaS)/sigmaT; % nullszórású
    
    %% analóg mintavételezés
    
    %tapasztalat: kihal mindig kivéve nagyon nagy hasadási hatkerre
    if m == 1 && (method ==0 || method == 1) % akkor fut le, ha mindhármat futtatom, vagy ha csak ezt
        quotientA = pF*nu + pS; % analitikus várható érték/sokszorozódás egy lépés után
        deviationA = (pF*nu^2 + pS - quotientA^2)/quotientA^2; %analitikus relatív szórásnégyzet
        for t = 1:tMax %idõléptetés
            if fP == 1
                if t == 1 %szabadúthossz
                    freePath(t) = - 1/sigmaT*log(rand);
                else
                    freePath(t) = freePath(t - 1) - 1/sigmaT*log(rand); % szabadúthossz, mennyivel haladjon elõre (1,1,1...)->(r,r,r...)
                end
            end
            for j = 1:nNeutron %neutronléptetés
                if rand < pA % abszorpció, most megszûnik, ha abszorbeál
                    w(j) =  0;
                elseif  rand < pF % hasadás, most 2.5 neutron keletkezik egy hasadásban
                    w(j) = nu*w(j);
                    f(t) = f(t) + 1; % hasadás + 1, hasadásszámláló
                end
            end
            nPop(t) = sum(w); % neutronpopuláció a t. idõben/lépés után
            var(t) = sum(w.^2)/sum(w)^2 - 1/nNeutron; % súlyok szórása adott lépés után
            Ew(t) = sum(w)/nNeutron; %súlyok várható értéke, átlagos súly adott lépés után
        end
        n = nPop(tMax); % neutronszám a szimuláció végén
        nFission = sum(f); % hasadások száma összesen
        deviationS = var(tMax);
        quotientS = nthroot(n/nNeutron,100); % szimuláció kvóciense, szimuláció által adott átlag szaporulat
    end
    
    %% implicit capture
    
    if m == 2 && (method == 0 || method == 2)
        quotientA = pF*nu + pS; % analitikus várható érték egy lépés után
        deviationA = (pF*nu^2 + pS - quotientA^2)/quotientA^2; %analitikus relatív szórásnégyzet
        for t = 1:tMax %idõléptetés
            if fP == 1
                if t == 1 %szabadúthossz
                    freePath(t) = - 1/sigmaT*log(rand);
                else
                    freePath(t) = freePath(t - 1) - 1/sigmaT*log(rand); % szabadúthossz, mennyivel haladjon elõre (1,1,1...)->(r,r,r...)
                end
            end
            for j = 1:nNeutron %neutronléptetés
                w(j) = (1 - pA)*w(j); %abszorpció valség, túlélés valségével lesúlyozom
                if  rand < pF/(pS + pF) % hasadás, most 2.5 neutron keletkezik egy hasadásban redukált eseménytérben
                    w(j) = nu*w(j);
                    f(t) = f(t) + 1; %hasadás számláló
                end
            end
            nPop(t) = sum(w); % neutronpopuláció/neutronteljesítmény a t. idõben
            var(t) = sum(w.^2)/sum(w)^2 - 1/nNeutron;
            Ew(t) = sum(w)/nNeutron; %súlyok várható értéke, átlagos súly
        end
        n = nPop(tMax);
        nFission = sum(f);
        deviationS = var(tMax);
        quotientS = nthroot(n/nNeutron,100); % szimuláció kvóciense, szimuláció által adott átlag szaporulat
    end
    
    %% nem analóg szimuláció
    
    if m == 3 && (method == 0 || method == 3)
        quotientA = pF*nu + pS; % analitikus várható érték egy lépés után
        deviationA = (pF*nu^2 + pS - quotientA^2)/quotientA^2; %analitikus relatív szórásnégyzet
        for t = 1:tMax %idõléptetés
            if fP == 1
                if t == 1 %szabadúthossz
                    freePath(t) = - 1/sigmaT*log(rand);
                else
                    freePath(t) = freePath(t - 1) - 1/sigmaT*log(rand); % szabadúthossz, mennyivel haladjon elõre (1,1,1...)->(r,r,r...)
                end
            end
            for j = 1:nNeutron %neutronléptetés
                w(j) = (1 - pA)*w(j); %túlélés valsége
                if rand  < q % hasadás, most 2.5 neutron keletkezik egy hasadásban
                    w(j) = 1/q*pF/(pF + pS)*nu*w(j);
                    f(t) = f(t) + 1;
                else
                    w(j) = w(j)*(1 - pF/(pF + pS))/(1 - q); % szóródik, súlykorrekció itt is kell
                end
            end
            nPop(t) = sum(w); % neutronpopuláció a t. idõben
            var(t) = sum(w.^2)/sum(w)^2 - 1/nNeutron;
            Ew(t) = sum(w)/nNeutron; %súlyok várható értéke, átlagos súly
        end
        n = nPop(tMax);
        nFission = sum(f);
        deviationS = var(tMax);
        quotientS = nthroot(n/nNeutron,tMax); % szimuláció kvóciense, szimuláció által adott átlag szaporulat 
    end

%% plotok
if method == 0 || (method > 0 && m == method)
figure(1)
hold on
if fP == 1
plot(freePath,nPop,'.-','MarkerSize',10,'LineWidth',1) % szabadúthossz függvényében
elseif fP == 0
plot(nPop,'.-','MarkerSize',10,'LineWidth',1) % egységnyi idõlépetetés függvényében
end
grid on
xlabel('Time')
ylabel('Neutron population')
title('Time dependence of neutron population')
end
%figure(2)
%plot(var,'b.-')
%grid on
%xlabel('Time')
%ylabel('Relative deviation')
%title('Time dependence of relative variance of weights')

%figure(3)
%plot(Ew,'g.-')
%grid on
%xlabel('idõ')
%ylabel('Average weight')
%title('Time dependence of average weight')
end
else
    disp('Rossz inputok!');
end
%legend('Simulated neutron pop.','Simulated neutron pop. 2','Simulated neutron pop. 3')
if fP == 0
fplot(@(t) nNeutron*quotientA^t,[0,tMax],'k');
end
if method == 0
legend('Analog','Implicit capture','Non analog','Analytic solution')
end
end