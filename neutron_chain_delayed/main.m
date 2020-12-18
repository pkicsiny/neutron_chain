function [ n, nFission, quotientA] = main( sigmaA,sigmaF,sigmaS,nNeutron,tMax)
    % nem anal�g szimul�ci� sz�r�sminimaliz�l�s
warning('off','all')
close all

nu = 2.5; % hasad�si sokszoroz�d�s
sigmaT = sigmaA + sigmaS + sigmaF;

pA=sigmaA/sigmaT;  %vals�gek [0,1]
pF=sigmaF/sigmaT;
pS=sigmaS/sigmaT;
for i = 1:200
    w = ones(1,nNeutron); % neutronok s�lyai
    f = zeros(1,tMax); % hasad�sok sz�ma
    nPop = zeros(1,tMax); % neutronpopul�ci� az id�ben
    var = zeros(1,tMax); % relat�v sz�r�sok vektora
    q(i) = 0 + 0.005*i;
    quotientA(i) = pF*nu + pS;
    %q(i) = 0.5+0.0005*i;
    for t = 1:tMax %id�l�ptet�s
        for j = 1:nNeutron %neutronl�ptet�s
            w(j) = (1 - pA)*w(j); %t�l�l�s vals�ge
            if rand  < q(i) % hasad�s, most 2.5 neutron keletkezik egy hasad�sban
                w(j) = 1/q(i)*pF/(pF + pS)*nu*w(j);
                f(t) = f(t) + 1;
            else
                w(j) = w(j)*(1 - pF/(pF + pS))/(1 - q(i)); % sz�r�dik, s�lykorrekci� itt is kell
            end
        end
        nPop(t) = sum(w); % neutronpopul�ci� a t. id�ben
        var(t) = sum(w.^2)/sum(w)^2 - 1/nNeutron; %s�lyok sz�r�sa
        Ew(t) = sum(w)/nNeutron; %s�lyok v�rhat� �rt�ke, �tlagos s�ly
    end
    n(i) = nPop(tMax); % neutronok sz�ma a szimul�ci� v�g�n
    nFission(i) = sum(f); % hasad�sok sz�ma
    qvar(i) = var(tMax); %szimul�ci� v�g�n a s�lyok sz�r�sa
    %varP(i) = sum(nPop.^2)/sum(nPop)^2 - 1/nNeutron; % neutronpopul�ci�
    quotientS(i) = nthroot(n(i)/nNeutron,tMax); % szimul�ci� kv�ciense, szimul�ci� �ltal adott �tlag szaporulat
    
end
figure(1)
plot(q,qvar,'b.-') %s�lyok sz�r�sn�gyzete
%plot(q,sqrt(qvar),'b.-') %s�lysz�r�s
grid on
xlabel('q')
ylabel('Relative deviation')
title('Relative deviation - q')
legend('Relative variance of the test')


figure(2)
hold on
plot(q,nFission,'k.-') % �sszes hasad�ssz�m egy szimul�ci�ban
plot(q,q*nNeutron*tMax,'r') %�sszes hasad�ssz�m analitikus �rt�ke
grid on
xlabel('q')
ylabel('Sum of fissions')
title('Sum of fissions - q')
hold off

figure(3)
hold on
plot(q,quotientS,'g.-') % szimul�ci� alapj�n sz�molt �tlagkv�ciens 
plot(q,quotientA,'r') % analitikus kv�ciens (v�rhat� �rt�ke egy l�p�s ut�n)
grid on
legend('Simulated quotient','Analytic quotient')
xlabel('q')
ylabel(' Expected population growth in one step (quotient)')
title('Quotient - q')
hold off 

figure(4)
hold on
plot(q,n,'m.-') % v�gs� neutronsz�m
plot(q,nNeutron*quotientA.^tMax,'r'); %elm�leti v�gs� neutronsz�m
grid on
legend('Simulated final population','Analytic final population')
xlabel('q')
ylabel('Final population')
title('Final popultion - q')
hold off


end

