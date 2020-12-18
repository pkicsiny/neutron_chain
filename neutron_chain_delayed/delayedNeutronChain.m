function [w,W,SumP,varq,lefut] = delayedNeutronChain( sigmaA,sigmaF,sigmaS,nNeutron,tMax,resolution,iLim)
    % paramétertorzításos szimuláció szórásminimalizálás életidõsorsolással
    %tMax - maximum idõ
    %resolution - max. idõt hányasával léptesse (tMax osztható legyen vele)
    %iLim = 0 - sima analóg
    %iLim = 1 - szakaszhatárra visszaállíás túllépés esetén
warning('off','all')
close all

if ~(tMax > 0 && resolution > 0 && sigmaA >= 0 && sigmaF >= 0 && sigmaS >= 0 && nNeutron > 0 && (iLim == 0 || iLim == 1 || iLim == 2))
error('Nem megfelelõ input paraméterek!');
exit(0);
end
%kérdés majd??????
betai=[21 142 128 257 75 27]*1e-5;
T12i=[55.7 22.7 6.2 2.3 0.615 0.23];
beta=sum(betai);
lambda=1/( sum(T12i.*betai)/beta );

nu = 2.5; % hasadási sokszorozódás
iMax = 300; % q(i), és q2(i) paraméter függésének felbontása
sigmaT = sigmaA + sigmaS + sigmaF;
pA=sigmaA/sigmaT;  %valségek [0,1]
pF=sigmaF/sigmaT;
pS=sigmaS/sigmaT;
if rem(tMax,resolution) == 0
interval = [resolution:resolution:tMax]; % idõintervallumok
else
    error('A max idõ legyen osztható a felbontással!');
    exit(0);
end
varq = zeros(length(interval),iMax); % relatív szórások mátrixa
SumP = zeros(iMax,length(interval));
lefut = 0;
bel = 0;
maxfel = 0;
for i = 1:iMax
    %minden q(i)-ra más:
    wT=zeros(length(interval),nNeutron);
    w = ones(1,nNeutron); % neutronok súlyai
    W = zeros(length(interval),nNeutron); % súlyok bankja
    if iLim == 2
        q(i) = (sigmaT - sigmaS)/sigmaT;
        q2(i) = 0.1 + 0.0016*i;
    else
        q(i) = 0.1 + 0.002*i;
    end
    for j = 1:nNeutron %neutronléptetés
        t = 0;
        k = 1; % intervallum léptetõ 
        while t < max(interval)
            w(j) = (1 - pA)*w(j); %túlélés valsége, implicit capture
            if rand  < q(i) % hasadás,prompt neutronok csak
                w(j) = pF/(pF + pS)/q(i)*w(j); 
                SumP(i,k) = SumP(i,k) + w(j); % leadott összteljesítmény (hasadások száma, súlyok összege) az adott intervallumban,
                % minden neutron által, mert a többi neutroné is ehhez fog
                % hozzáadódni (neutronpopuláció)
                W(k,j) = W(k,j) + w(j); % leadott teljesítmény neutrononként, intervallumonként eltárolva: sum(W(k,:)) = SumP(k) teljesül
                w(j) = w(j)*nu;
                dt=-1/lambda*log(rand); % életidõ sorsolás (sigmaT vagy lambda???????)
                if iLim == 2
                dt=-1/q2(i)*log(rand); % életidõ sorsolás torzítottan
                end
                t = t + dt;
                if t > interval(k) %&& ter < interval(k)  % ha az idõ túllépi az adott intervallumot, akkor a következõbe lép
                    if iLim == 1 % szakaszhatárra visszaállítja az idõt
                        t = interval(k);
                        dt = -1/lambda*log(rand); % új életidõ sorsolás dt2
                        t = t + dt;
                    elseif iLim == 2 % visszaállít + visszatorzít
                        w1 = exp(-(interval(k) - (t - dt))*(lambda - q2(i))); %(t - dt) helyett t volt
                        w(j) = w(j)*w1; %súlykorrekció
                        t = interval(k); % szakaszhatárra állítás
                        wT(k,j) = w(j); % elmentjük
                        maxfel= maxfel + 1;
                        if t < max(interval)  %!!!!!! <= nem kell csak <, mert akkor az utolsónál is belép, de ott nem kell
                         %az utsó intervallumnál nem kell   
                        dt = -1/q2(i)*log(rand);
                        t = t + dt; % ha dt2 > resolution (5), k kettõt lép!!
                        w1 = lambda/q2(i)*exp(-dt*(lambda - q2(i)));
                        w(j) = w(j)*w1;
                        lefut = lefut + 1;
                        end
                    end
                    if k + 1 + floor(dt/resolution) <= length(interval)%ha nem ugrik ki
                    k = k + 1 + floor(dt/resolution); % egynél többet lép itt esetenként! nem 5tel hanem a lépéssel kell leosztani
                    else % (10,20),30,40 re is ugorhat, ilyenkor 20ra kell állítani, különben k index túllép
                        k = length(interval);
                    end
                    elseif iLim == 2 && t < interval(k) % ha nem lépi túl az intervallumhatárt, akkor is kell korrekció
                   
                    w0 = lambda/q2(i)*exp(-dt*(lambda - q2(i)));
                    w(j) = w(j)*w0;
                    bel = bel + 1;
                end
            else  % szóródik, súlykorrekció itt is kell
                w(j) = w(j)*(pS/(pF + pS))/(1 - q(i));
            end
        end
    end
    for k = 1:length(interval)
        Ws=W(k,W(k,:)>0); hossz=length(Ws); %1 volt mindig az 1. k helyett:(
        varq(k,i) = sqrt(sum(Ws.^2)/(sum(Ws)^2) - 1/nNeutron); % 0-k nincsenek benne, rövidebb
    % varq(k,i) = sqrt(sum(W(k,:).^2)/(sum(W(k,:))^2) - 1/nNeutron); %súlyok szórásnégyzete a k. intervallumban q(i) paraméter esetén
    end
end

figure(1)
hold on
for k = 1:length(interval)
    if iLim == 2
plot(q2,varq(k,:),'.-') %súlyok szórása
xlabel('q2 parameter')
    else
        plot(q,varq(k,:),'.-') %súlyok szórása
        xlabel('q parameter')
    end
end
grid on
ylabel('Relative deviation')
title('Relative deviation - q2')

figure(2)
hold on
if iLim == 2
    surf(q2,interval,varq)
     xlabel('q2')
else
    surf(q,interval,varq)
     xlabel('q')
end
grid on
% ahol nincs ott NaN a szórás, a  szélsõséges paramétertorzítás miatt, azaz kihalt
%pontokat kell figyelni, interval nem a 0-tól indul, az elsõ pont az 1.
%intervallumban a szórás
 ylabel('Time interval')
 zlabel('Relative deviation')
 title('parameter dependence of the relative variance in each interval')
hold off

figure(3)
hold on
if iLim == 2
    surf(q2,interval,log(SumP'))
   % surf(q2,interval,log10(SumP'))
    xlabel('q2')
else
    surf(q,interval,SumP')
    %surf(q,interval,log10(SumP'))
    xlabel('q')
end
grid on
 ylabel('Time interval')
 zlabel('Reactor power')
 title('parameter dependence of the reactor powers in each interval')
hold off
end