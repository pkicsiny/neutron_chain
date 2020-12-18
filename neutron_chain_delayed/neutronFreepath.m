function [ n, nFission,quotientA,quotientS,deviationA,deviationS] = neutronFreePath( sigmaA,sigmaF,sigmaS,nNeutron,tMax,fP,method )
%tMax - maximum id�l�ptet�s
%fP = 0 - nincs szabad�thossz sorsol�s 
%fP = 1 - van
% method = 0 - mindh�rom m�dszer lefut
% method = 1/2/3 - 1./2./3. m�dszer fut csak le

warning('off','all')
close all
if (fP == 0 || fP == 1) && (method == 0 || method == 1 || method == 2 || method == 3) && nNeutron > 0 && tMax > 0 && sigmaA >= 0 && sigmaF >= 0 && sigmaS >= 0 
for m = 1:3
    nu = 2.5; % hasad�si sokszoroz�d�s
    sigmaT = sigmaA + sigmaS + sigmaF; % teljes hat�skeresztmetszet
    w = ones(1,nNeutron); % neutronok s�lyai (1,1,1,....) kezdetben
    f = zeros(1,tMax); % hasad�sok sz�ma adott l�p�sben
    nPop = zeros(1,tMax); % neutronpopul�ci�/teljes�tm�ny az id�ben
    var = zeros(1,tMax); % s�lyok relat�v sz�r�s�nak vektora
    pA=sigmaA/sigmaT;  %vals�gek [0,1]
    pF=sigmaF/sigmaT;
    pS=sigmaS/sigmaT;
    freePath = zeros(1,tMax); % szabad�thossz (id� most????)
    %q = (nu*pF)/(pS + nu*pF); %ezzel a s�lysz�r�s minimaliz�l�dik, n�pess�gszab�lyoz�s
    q = (sigmaT-sigmaS)/sigmaT; % nullsz�r�s�
    
    %% anal�g mintav�telez�s
    
    %tapasztalat: kihal mindig kiv�ve nagyon nagy hasad�si hatkerre
    if m == 1 && (method ==0 || method == 1) % akkor fut le, ha mindh�rmat futtatom, vagy ha csak ezt
        quotientA = pF*nu + pS; % analitikus v�rhat� �rt�k/sokszoroz�d�s egy l�p�s ut�n
        deviationA = (pF*nu^2 + pS - quotientA^2)/quotientA^2; %analitikus relat�v sz�r�sn�gyzet
        for t = 1:tMax %id�l�ptet�s
            if fP == 1
                if t == 1 %szabad�thossz
                    freePath(t) = - 1/sigmaT*log(rand);
                else
                    freePath(t) = freePath(t - 1) - 1/sigmaT*log(rand); % szabad�thossz, mennyivel haladjon el�re (1,1,1...)->(r,r,r...)
                end
            end
            for j = 1:nNeutron %neutronl�ptet�s
                if rand < pA % abszorpci�, most megsz�nik, ha abszorbe�l
                    w(j) =  0;
                elseif  rand < pF % hasad�s, most 2.5 neutron keletkezik egy hasad�sban
                    w(j) = nu*w(j);
                    f(t) = f(t) + 1; % hasad�s + 1, hasad�ssz�ml�l�
                end
            end
            nPop(t) = sum(w); % neutronpopul�ci� a t. id�ben/l�p�s ut�n
            var(t) = sum(w.^2)/sum(w)^2 - 1/nNeutron; % s�lyok sz�r�sa adott l�p�s ut�n
            Ew(t) = sum(w)/nNeutron; %s�lyok v�rhat� �rt�ke, �tlagos s�ly adott l�p�s ut�n
        end
        n = nPop(tMax); % neutronsz�m a szimul�ci� v�g�n
        nFission = sum(f); % hasad�sok sz�ma �sszesen
        deviationS = var(tMax);
        quotientS = nthroot(n/nNeutron,100); % szimul�ci� kv�ciense, szimul�ci� �ltal adott �tlag szaporulat
    end
    
    %% implicit capture
    
    if m == 2 && (method == 0 || method == 2)
        quotientA = pF*nu + pS; % analitikus v�rhat� �rt�k egy l�p�s ut�n
        deviationA = (pF*nu^2 + pS - quotientA^2)/quotientA^2; %analitikus relat�v sz�r�sn�gyzet
        for t = 1:tMax %id�l�ptet�s
            if fP == 1
                if t == 1 %szabad�thossz
                    freePath(t) = - 1/sigmaT*log(rand);
                else
                    freePath(t) = freePath(t - 1) - 1/sigmaT*log(rand); % szabad�thossz, mennyivel haladjon el�re (1,1,1...)->(r,r,r...)
                end
            end
            for j = 1:nNeutron %neutronl�ptet�s
                w(j) = (1 - pA)*w(j); %abszorpci� vals�g, t�l�l�s vals�g�vel les�lyozom
                if  rand < pF/(pS + pF) % hasad�s, most 2.5 neutron keletkezik egy hasad�sban reduk�lt esem�nyt�rben
                    w(j) = nu*w(j);
                    f(t) = f(t) + 1; %hasad�s sz�ml�l�
                end
            end
            nPop(t) = sum(w); % neutronpopul�ci�/neutronteljes�tm�ny a t. id�ben
            var(t) = sum(w.^2)/sum(w)^2 - 1/nNeutron;
            Ew(t) = sum(w)/nNeutron; %s�lyok v�rhat� �rt�ke, �tlagos s�ly
        end
        n = nPop(tMax);
        nFission = sum(f);
        deviationS = var(tMax);
        quotientS = nthroot(n/nNeutron,100); % szimul�ci� kv�ciense, szimul�ci� �ltal adott �tlag szaporulat
    end
    
    %% nem anal�g szimul�ci�
    
    if m == 3 && (method == 0 || method == 3)
        quotientA = pF*nu + pS; % analitikus v�rhat� �rt�k egy l�p�s ut�n
        deviationA = (pF*nu^2 + pS - quotientA^2)/quotientA^2; %analitikus relat�v sz�r�sn�gyzet
        for t = 1:tMax %id�l�ptet�s
            if fP == 1
                if t == 1 %szabad�thossz
                    freePath(t) = - 1/sigmaT*log(rand);
                else
                    freePath(t) = freePath(t - 1) - 1/sigmaT*log(rand); % szabad�thossz, mennyivel haladjon el�re (1,1,1...)->(r,r,r...)
                end
            end
            for j = 1:nNeutron %neutronl�ptet�s
                w(j) = (1 - pA)*w(j); %t�l�l�s vals�ge
                if rand  < q % hasad�s, most 2.5 neutron keletkezik egy hasad�sban
                    w(j) = 1/q*pF/(pF + pS)*nu*w(j);
                    f(t) = f(t) + 1;
                else
                    w(j) = w(j)*(1 - pF/(pF + pS))/(1 - q); % sz�r�dik, s�lykorrekci� itt is kell
                end
            end
            nPop(t) = sum(w); % neutronpopul�ci� a t. id�ben
            var(t) = sum(w.^2)/sum(w)^2 - 1/nNeutron;
            Ew(t) = sum(w)/nNeutron; %s�lyok v�rhat� �rt�ke, �tlagos s�ly
        end
        n = nPop(tMax);
        nFission = sum(f);
        deviationS = var(tMax);
        quotientS = nthroot(n/nNeutron,tMax); % szimul�ci� kv�ciense, szimul�ci� �ltal adott �tlag szaporulat 
    end

%% plotok
if method == 0 || (method > 0 && m == method)
figure(1)
hold on
if fP == 1
plot(freePath,nPop,'.-','MarkerSize',10,'LineWidth',1) % szabad�thossz f�ggv�ny�ben
elseif fP == 0
plot(nPop,'.-','MarkerSize',10,'LineWidth',1) % egys�gnyi id�l�petet�s f�ggv�ny�ben
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
%xlabel('id�')
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