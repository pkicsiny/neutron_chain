function [w,W,wT,SumP,varq,varq2,varq3,lambda,beta,ujlambda,qopt,betaopt,lambdaopt,quotient] = betaOpt(pF,pS,nNeutron,tMax,resolution)
% paramétertorzításos szimuláció szórásminimalizálás életidõsorsolással
%tMax - maximum idõ
%resolution - max. idõt hányasával léptesse (tMax osztható legyen vele)
%iLim = 0 - sima analóg
%iLim = 1 - szakaszhatárra visszaállíás túllépés esetén
warning('off','all')
close all
pA = 1 - pF - pS;
if ~(tMax > 0 && resolution > 0 && pA >= 0 && pF >= 0 && pS >= 0 && nNeutron > 0 && pF <= 1 && pS <= 1)
    error('Nem megfelelõ input paraméterek!');
    exit(0);
end
%kérdés majd??????
betai=[21 142 128 257 75 27]*1e-5;
T12i=[55.7 22.7 6.2 2.3 0.615 0.23];
beta=sum(betai);
lambda=1/( sum(T12i.*betai)/beta );

nu = 2.5; % hasadási sokszorozódás
iMax = 40; % q(i), és q2(i) paraméter függésének felbontása
bMax = 40; % új béta paraméter függésének felbontása
%ujbeta = 0.3;
if rem(tMax,resolution) == 0
    interval = resolution:resolution:tMax; % idõintervallumok
else %hibaüzenet
    error('A max idõ legyen osztható a felbontással!');
    exit(0);
end
varq = zeros(bMax,iMax); % relatív szórások mátrixa
SumP = zeros(iMax,bMax);
lefut = 0;
bel = 0;
maxfel = 0;

qopt = 1 - pS;
betaopt = 1 - pF*nu*(1 - beta)/qopt;
lambdaopt = lambda*pF*nu*beta/(qopt*betaopt);
quotient = pS + (1 - beta)*pF*nu;

for i = 1:iMax
    q(i) = qopt;
    ujlambda(i) = 0.0805 + 0.005*i;
    for b = 1:bMax
        %minden q(i)-ra más:
        wT=zeros(length(interval),nNeutron);
        w = ones(1,nNeutron); % neutronok súlyai
        W = zeros(length(interval),nNeutron); % súlyok bankja
        ujbeta(b) =  0.0032 + 0.0001*b;
        for j = 1:nNeutron %neutronléptetés
            t = 0;
            %t = -1/lambda*log(rand); % kezdeti idõt sorsolunk minden neutronnak
            k = 1; % intervallum léptetõ
            while t < max(interval)
                w(j) = (1 - pA)*w(j); %túlélés valsége, implicit capture
                if rand  < q(i) % hasadás,prompt neutronok csak
                    w(j) = pF/(pF + pS)/q(i)*w(j);
                    %most SumP = fiplus
                    SumP(i,b) = SumP(i,b) + w(j); % leadott összteljesítmény (hasadások száma, súlyok összege)(i,k) volt, most fiplus
                    W(k,j) = W(k,j) + w(j); % leadott teljesítmény neutrononként, intervallumonként eltárolva: sum(W(k,:)) = SumP(k) teljesül
                    w(j) = w(j)*nu;
                    
                    if rand < ujbeta(b) % itt keletkezik késsõneutron: ez csak idõben ugrik
                        w(j)=w(j)*beta/ujbeta(b); % uj béta korrekció
                        dt=-1/ujlambda(i)*log(rand); % életidõ sorsolás torzítottan
                        t = t + dt;
                        if t > interval(k) % ha az idõ túllépi az adott intervallumot, akkor a következõbe lép
                            w1 = exp(-(interval(k) - (t - dt))*(lambda - ujlambda(i))); %(t - dt) helyett t volt
                            w(j) = w(j)*w1; %súlykorrekció
                            t = interval(k); % szakaszhatárra állítás
                            wT(k,j) = w(j); % az adott intervallum határON a neutronlánc súlya
                            maxfel= maxfel + 1;
                            if t < max(interval)  %!!!!!! <= nem kell csak <, mert akkor az utolsónál is belép, de ott nem kell
                                %az utsó intervallumnál nem kell
                                dt = -1/ujlambda(i)*log(rand);
                                t = t + dt; % ha dt2 > resolution (5), k kettõt lép!!
                                w1 = lambda/ujlambda(i)*exp(-dt*(lambda - ujlambda(i)));
                                w(j) = w(j)*w1;
                                lefut = lefut + 1;
                            end
                            if k + 1 + floor(dt/resolution) <= length(interval)%ha nem ugrik ki
                                k = k + 1 + floor(dt/resolution); % egynél többet lép itt esetenként! nem 5tel hanem a lépéssel kell leosztani
                            else % (10,20),30,40 re is ugorhat, ilyenkor 20ra kell állítani, különben k index túllép
                                k = length(interval);
                            end
                        elseif  t < interval(k) % ha nem lépi túl az intervallumhatárt, akkor is kell korrekció
                            
                            w0 = lambda/ujlambda(i)*exp(-dt*(lambda - ujlambda(i)));
                            w(j) = w(j)*w0;
                            bel = bel + 1;
                        end
                    else
                        w(j) = w(j)*(1 - beta)/(1 - ujbeta(b)); % béta korr.
                    end
                    
                else  % szóródik, súlykorrekció itt is kell
                    w(j) = w(j)*(pS/(pF + pS))/(1 - q(i));
                end
            end
        end
        %for k = 1:length(interval) %szórások, most csak a végsõ intervallumot fogja nézni béta' és lambda' függvényében
        Ws=W(length(interval),W(length(interval),:)>0); %hossz=length(Ws); %1 volt mindig az 1. k helyett:(
        varq(b,i) = sqrt(sum(Ws.^2)/(sum(Ws)^2) - 1/nNeutron); % 0-k nincsenek benne, rövidebb
        varq2(b,i) = var(wT(length(interval),wT(length(interval),:)>0)); % sima szórás a végén
        Ws3=wT(length(interval),wT(length(interval),:)>0); %hossz=length(Ws3); %relatív szórás a végén
        varq3(b,i) = sqrt(sum(Ws3.^2)/(sum(Ws3)^2) - 1/nNeutron);
        % end
        
    end
end

figure(1)
hold on
surf(ujlambda,ujbeta,varq)
xlabel('\lambda''','fontsize',18)
grid on
% ahol nincs ott NaN a szórás, a  szélsõséges paramétertorzítás miatt, azaz kihalt
%pontokat kell figyelni, interval nem a 0-tól indul, az elsõ pont az 1.
%intervallumban a szórás
ylabel('\beta''','fontsize',18)
zlabel('RSD','fontsize',18)
title({'Biasing parameter dependence of the RSD','of power output in final interval'},'fontsize',16)
hold off

figure(2)
hold on
surf(ujlambda,ujbeta,varq3)
xlabel('\lambda''','fontsize',18)
grid on
%intervallumok határán lévõ súlyok szórása
ylabel('\beta''','fontsize',18)
zlabel('RSD','fontsize',18)
title({'RSD of weight distribution','at the end of the final interval'},'fontsize',16)
hold off

figure(3)
hold on
%surf(ujlambda,ujbeta,SumP') %fiplus nem kell mert sumP már az egészet
%tartalmazza
surf(ujlambda,ujbeta,SumP')
set(gca,'zscale','log');
xlabel('\lambda''','fontsize',18)
grid on
ylabel('\beta''','fontsize',18)
zlabel('Reactor power output','fontsize',18)
title({'Parameter dependence of the','total power output'},'fontsize',16)
hold off

figure(4)
hold on
surf(ujlambda,ujbeta,varq2)
set(gca,'zscale','log');
xlabel('\lambda''','fontsize',18)
grid on
%intervallumok határán lévõ súlyok szórása
ylabel('\beta''','fontsize',18)
zlabel('Variance','fontsize',18)
title({'Variance of weight distribution','at the end of the final interval'},'fontsize',16)
hold off
end
