function [w,W,wT,SumP,varq,varq2,beta,fiplus,q2] = beta( sigmaA,sigmaF,sigmaS,nNeutron,tMax,resolution)
% paramétertorzításos szimuláció szórásminimalizálás életidõsorsolással
%tMax - maximum idõ
%resolution - max. idõt hányasával léptesse (tMax osztható legyen vele)
%iLim = 0 - sima analóg
%iLim = 1 - szakaszhatárra visszaállíás túllépés esetén
warning('off','all')
close all

if ~(tMax > 0 && resolution > 0 && sigmaA >= 0 && sigmaF >= 0 && sigmaS >= 0 && nNeutron > 0)
    error('Nem megfelelõ input paraméterek!');
    exit(0);
end
%kérdés majd??????
betai=[21 142 128 257 75 27]*1e-5;
T12i=[55.7 22.7 6.2 2.3 0.615 0.23];
beta=sum(betai);
lambda=1/( sum(T12i.*betai)/beta );

nu = 2.5; % hasadási sokszorozódás
iMax = 100; % q(i), és q2(i) paraméter függésének felbontása
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
    q(i) = (sigmaT - sigmaS)/sigmaT;
    q2(i) = 0.08 + 0.0008725*i;
    for j = 1:nNeutron %neutronléptetés
        t = 0;
        %t = -1/lambda*log(rand); % kezdeti idõt sorsolunk minden neutronnak
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
                %if rand < beta % ha késõneutron akkor idõben ugrik
                    dt=-1/q2(i)*log(rand); % életidõ sorsolás torzítottan
                    t = t + dt;
                    if t > interval(k) % ha az idõ túllépi az adott intervallumot, akkor a következõbe lép
                        w1 = exp(-(interval(k) - (t - dt))*(lambda - q2(i))); %(t - dt) helyett t volt
                        w(j) = w(j)*w1; %súlykorrekció
                        t = interval(k); % szakaszhatárra állítás
                        wT(k,j) = w(j); % az adott intervallum határON a neutronlánc súlya
                        maxfel= maxfel + 1;
                        if t < max(interval)  %!!!!!! <= nem kell csak <, mert akkor az utolsónál is belép, de ott nem kell
                            %az utsó intervallumnál nem kell
                            dt = -1/q2(i)*log(rand);
                            t = t + dt; % ha dt2 > resolution (5), k kettõt lép!!
                            w1 = lambda/q2(i)*exp(-dt*(lambda - q2(i)));
                            w(j) = w(j)*w1;
                            lefut = lefut + 1;
                        end
                        if k + 1 + floor(dt/resolution) <= length(interval)%ha nem ugrik ki
                            k = k + 1 + floor(dt/resolution); % egynél többet lép itt esetenként! nem 5tel hanem a lépéssel kell leosztani
                        else % (10,20),30,40 re is ugorhat, ilyenkor 20ra kell állítani, különben k index túllép
                            k = length(interval);
                        end
                    elseif  t < interval(k) % ha nem lépi túl az intervallumhatárt, akkor is kell korrekció
                        
                        w0 = lambda/q2(i)*exp(-dt*(lambda - q2(i)));
                        w(j) = w(j)*w0;
                        bel = bel + 1;
                    end
                %end
            else  % szóródik, súlykorrekció itt is kell
                w(j) = w(j)*(pS/(pF + pS))/(1 - q(i));
            end
        end
    end
    for k = 1:length(interval)
        Ws=W(k,W(k,:)>0); hossz=length(Ws); %1 volt mindig az 1. k helyett:(
        varq(k,i) = sum(Ws.^2)/(sum(Ws)^2) - 1/nNeutron; % 0-k nincsenek benne, rövidebb
        varq2(k,i) = var(wT(k,wT(k,:)>0)); % sima szórás a határokon
        Ws3=wT(k,wT(k,:)>0); hossz=length(Ws3); %határon súly rel szórás
        varq3(k,i) = sum(Ws3.^2)/(sum(Ws3)^2) - 1/nNeutron;
    end
end

fiplus = sum(SumP'); % össz leadott teljesítmény

figure(1)
hold on
surf(q2,interval,varq)
%surf(q2,interval,reference)
xlabel('\lambda''','fontsize',18)
grid on
% ahol nincs ott NaN a szórás, a  szélsõséges paramétertorzítás miatt, azaz kihalt
%pontokat kell figyelni, interval nem a 0-tól indul, az elsõ pont az 1.
%intervallumban a szórás
ylabel('Time interval','fontsize',18)
zlabel('RSD','fontsize',18)
title({'Biasing parameter dependence of the RSD','of power output in each interval'},'fontsize',16)
hold off

figure(2)
hold on
surf(q2,interval,varq2)
xlabel('\lambda''','fontsize',18)
grid on
%intervallumok határán lévõ súlyok szórása
ylabel('Time interval','fontsize',18)
zlabel('RSD','fontsize',18)
title({'Biasing parameter dependence of the RSD','of weight distribution at interval limits'},'fontsize',16)
hold off

figure(3)
hold on
surf(q2,interval,SumP')
 %surf(q2,interval,log10(SumP'))
xlabel('\lambda''','fontsize',18)
grid on
ylabel('Time interval','fontsize',18)
zlabel('Reactor power output','fontsize',18)
title({'Parameter dependence of the reactor',' power output in each interval'},'fontsize',16)
hold off

figure(4)
hold on
for k = 1:length(interval)
    plot(q2,varq(k,:),'.-') %súlyok szórása az intervallumban
    xlabel('\lambda'' parameter','fontsize',18)
end
grid on
ylabel('RSD','fontsize',18)
title({'Biasing parameter dependence of the RSD','of power output in each interval'},'fontsize',16)
legend('t<t1','t1<t<t2','t2<t<t3','t3<t<t4')

figure(5)
hold on
for k = 1:length(interval)
    plot(q2,varq2(k,:),'.-') %súlyok szórása a határokon
    xlabel('\lambda'' parameter','fontsize',18)
end
grid on
ylabel('RSD','fontsize',18)
title({'Biasing parameter dependence of the RSD','of weight distribution at interval limits'},'fontsize',16)
legend('t=t1','t=t2','t=t3','t=t4')

figure(6)
hold on
plot(q2,fiplus)
xlabel('\lambda''','fontsize',18)
grid on
ylabel('Total power output','fontsize',18)
title({'Parameter dependence of the total',' power output'},'fontsize',16)
hold off
end