function [w,W,wT,SumP,varq,lefut,fiplus] = phibeta( sigmaA,sigmaF,sigmaS,nNeutron,tMax,resolution)
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
betai=[21 142 128 257 75 27]*1e-5;
T12i=[55.7 22.7 6.2 2.3 0.615 0.23];
beta=sum(betai);
lambda=1/( sum(T12i.*betai)/beta );

nu = 2.5; % hasadási sokszorozódás
iMax = 200; % q(i), és q2(i) paraméter függésének felbontása
sigmaT = sigmaA + sigmaS + sigmaF;
pA=sigmaA/sigmaT;  %valségek [0,1]
pF=sigmaF/sigmaT;
pS=sigmaS/sigmaT;
if rem(tMax,resolution) == 0
    interval = resolution:resolution:tMax; % idõintervallumok
else
    error('A max idõ legyen osztható a felbontással!');
    exit(0);
end

fiplus = zeros(length(interval) + 1,iMax); % jövõbeli teljesítmények adott t-tõl kezdve: elsõ sor a teljes, utolsó sor 0
fiplus(length(interval) + 1,:) = 0; %utolsó 0 pl. 20 a max idõ és 20-tól indítom (csak a szépség kedvéért)
init = [0,interval];%(1:end-1)]; [0,5,10,15,20,..]

for z = 1:length(interval) % jövõbeli teljesítmény, szûzen indítom mindig
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
    q2(i) = 0.083625 + 0.00083625*i;
    for j = 1:nNeutron %neutronléptetés
        t = init(z); %% most mindig késõbbrõl indul  |__|__|__|__|
        k = z; % intervallum léptetõ %%1, most mindig késõbbrõl indul
        while t < max(interval)
            w(j) = (1 - pA)*w(j); %túlélés valsége, implicit capture
            if rand  < q(i) % hasadás,prompt neutronok csak
                w(j) = pF/(pF + pS)/q(i)*w(j);
                SumP(i,k) = SumP(i,k) + w(j); % leadott összteljesítmény (hasadások száma, súlyok összege) az adott intervallumban,
                % minden neutron által, mert a többi neutroné is ehhez fog
                % hozzáadódni (neutronpopuláció)
                W(k,j) = W(k,j) + w(j); % leadott teljesítmény neutrononként, intervallumonként eltárolva: sum(W(k,:)) = SumP(k) teljesül
                w(j) = w(j)*nu;
                
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
                
            else  % szóródik, súlykorrekció itt is kell
                w(j) = w(j)*(pS/(pF + pS))/(1 - q(i));
            end
        end
    end
    for k = 1:length(interval)
        Ws=W(k,W(k,:)>0); %hossz=length(Ws); %1 volt mindig az 1. k helyett:(
        varq(k,i) = sum(Ws.^2)/(sum(Ws)^2) - 1/nNeutron; % 0-k nincsenek benne, rövidebb
        varq2(k,i) = var(wT(k,wT(k,:)>0)); % sima szórás a határokon
        Ws3=wT(k,wT(k,:)>0); %hossz=length(Ws3); %határon súly rel szórás
        varq3(k,i) = sum(Ws3.^2)/(sum(Ws3)^2) - 1/nNeutron;
    end
end
fiplus(z,:) = sum(SumP'); % össz leadott teljesítmény
end

figure(1)
hold on
surf(q2,init,fiplus)
xlabel('\lambda''','fontsize',18)
grid on
%intervallumok határán lévõ súlyok szórása
ylabel('Time interval','fontsize',18)
zlabel('Total power output','fontsize',18)
title({'Biasing parameter and time dependence of the total','power output'},'fontsize',16)
hold off
end
