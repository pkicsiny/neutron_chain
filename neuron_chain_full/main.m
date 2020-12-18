function [ n, W,SumP, quotientA,w] = main( sigmaA,sigmaF,sigmaS,nNeutron,tMax)
    % nem anal�g szimul�ci� sz�r�sminimaliz�l�s
warning('off','all')
close all

nu = 2.5; % hasad�si sokszoroz�d�s
sigmaT = sigmaA + sigmaS + sigmaF;

pA=sigmaA/sigmaT;  %vals�gek [0,1]
pF=sigmaF/sigmaT;
pS=sigmaS/sigmaT;

SumP = zeros(1,200);
W= zeros(1,nNeutron);
for i = 1:200 % i a q-kra megy
    W= zeros(1,nNeutron); % hasad�sok sz�ma (fi)
    w = ones(1,nNeutron); % neutronok s�lyai (popul�ci�kontroll)
    nPop = zeros(1,tMax); % neutronpopul�ci� az id�ben
    var = zeros(1,tMax); % relat�v sz�r�sok vektora (pop kontroll)
    varW = zeros(1,tMax); % hasad�sok sz�r�sa (fi)
    %q(i) = 0 + 0.005*i; %teljes t�rk�p
    quotientA(i) = pF*nu + pS;
     %q(i) = 0.5+0.0005*i; % pop kontroll k�zeli
     q(i) = 0.5+0.0005*i; % fi k�zeli
    for t = 1:tMax %id�l�ptet�s
        for j = 1:nNeutron %neutronl�ptet�s
            w(j) = (1 - pA)*w(j); %t�l�l�s vals�ge
            if rand  < q(i) % hasad�s, most 2.5 neutron keletkezik egy hasad�sban
                w(j) = 1/q(i)*pF/(pF + pS)*w(j);
                SumP(i) = SumP(i) + w(j);
                W(j) = W(j) + w(j);
                w(j) = nu*w(j);
            else
                w(j) = w(j)*(1 - pF/(pF + pS))/(1 - q(i)); % sz�r�dik, s�lykorrekci� itt is kell
            end
        end
        nPop(t) = sum(w); % neutronpopul�ci� a t. id�ben
        var(t) = sum(w.^2)/sum(w)^2 - 1/nNeutron; %s�lyok sz�r�sa
        varW(t) = sum(W.^2)/sum(W)^2 - 1/nNeutron; %s�lyok sz�r�sa
        Ew(t) = sum(w)/nNeutron; %s�lyok v�rhat� �rt�ke, �tlagos s�ly
    end
    n(i) = nPop(tMax); % neutronok sz�ma a szimul�ci� v�g�n
    qvar(i) = var(tMax); %szimul�ci� v�g�n a s�lyok sz�r�sa
    qvarW(i) = varW(tMax); %szimul�ci� v�g�n a s�lyok sz�r�sa
    quotientS(i) = nthroot(n(i)/nNeutron,tMax); % szimul�ci� kv�ciense, szimul�ci� �ltal adott �tlag szaporulat
    
end
figure(1)
plot(q,qvar,'b.-') %s�lyok sz�r�sn�gyzete
%plot(q,sqrt(qvar),'b.-') %s�lysz�r�s
grid on
xlabel('q','fontsize',18)
ylabel('Relative variance (RSD^2)','fontsize',18)
title('RSD^2 - q','fontsize',18)
legend('Relative variance of the neutron weights','fontsize',18)

figure(2)
%plot(q,qvarW,'r.-') %s�lyok sz�r�sn�gyzete
plot(q,sqrt(qvarW),'r.-') %s�lysz�r�s
grid on
xlabel('q','fontsize',18)
ylabel('Relative variance (RSD^2)','fontsize',18)
title('RSD^2 - q','fontsize',18)
legend('Relative variance of the power output','fontsize',18)

figure(3)
hold on
plot(q,quotientS,'g.-') % szimul�ci� alapj�n sz�molt �tlagkv�ciens 
plot(q,quotientA,'r') % analitikus kv�ciens (v�rhat� �rt�ke egy l�p�s ut�n)
grid on
legend('Simulated quotient','Analytic quotient','fontsize',18)
xlabel('q','fontsize',18)
ylabel('Average population growth in one step','fontsize',14)
title('Quotient - q','fontsize',18)
hold off 

figure(4)
hold on
plot(q,n,'m.-') % v�gs� neutronsz�m
plot(q,nNeutron*quotientA.^tMax,'r'); %elm�leti v�gs� neutronsz�m
grid on
legend('Simulated final population','Analytic final population','fontsize',18)
xlabel('q','fontsize',18)
ylabel('Final population','fontsize',18)
title('Final popultion - q','fontsize',18)
hold off


end

