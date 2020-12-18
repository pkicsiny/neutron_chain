warning('off','all')
close all
pF = 0.298;
pS = 0.256;
nNeutron = 1000;
tMax = 100;
resolution = 1;
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
if rem(tMax,resolution) == 0
    interval = resolution:resolution:tMax; % idõintervallumok
else %hibaüzenet
    error('A max idõ legyen osztható a felbontással!');
    exit(0);
end
qopt = 1 - pS;
betaopt = 1 - pF*nu*(1 - beta)/qopt;
lambdaopt = lambda*pF*nu*beta/(qopt*betaopt);

quotient = pS + (1 - beta)*pF*nu;
szorzo = zeros(length(interval) + 1,1);
integral = zeros(length(interval),1);
fiplus = zeros(length(interval) + 1,1); % jövõbeli teljesítmények adott t-tõl kezdve: elsõ sor a teljes, utolsó sor 0
fiplus(length(interval) + 1) = 0; %utolsó 0 pl. 20 a max idõ és 20-tól indítom (csak a szépség kedvéért)
init = [0,interval];%(1:end-1)]; [0,5,10,15,20,..]
%SumPossz = zeros(length(interval),1);
s = 0;
for z = 1:length(interval)
varq = zeros(length(interval),1); % relatív szórások mátrixa
SumP = zeros(length(interval),nNeutron);
SumPossz = 0;
Egyszer = 0;
lefut = 0;
        wT=zeros(length(interval),nNeutron);
        w = ones(1,nNeutron); % neutronok súlyai
        W = zeros(length(interval),nNeutron); % súlyok bankja
        for j = 1:nNeutron %neutronléptetés
            t = init(z);
           % t = 0;
            s2 = 0;
            k = z; % intervallum léptetõ
           %k = 1;
            while t < max(interval)
                w(j) = (1 - pA)*w(j); %túlélés valsége, implicit capture
                if rand  < qopt % hasadás,prompt neutronok csak
                    w(j) = pF/(pF + pS)/qopt*w(j);
                    SumP(k,j) = SumP(k,j) + w(j);
                    SumPossz = SumPossz + w(j); % leadott összteljesítmény (hasadások száma, súlyok összege)(i,k) volt
                    W(k,j) = W(k,j) + w(j); % leadott teljesítmény neutrononként, intervallumonként eltárolva: sum(W(k,:)) = SumP(k) teljesül
                    w(j) = w(j)*nu;
                    r = rand;
                    if r < betaopt % itt keletkezik késsõneutron: ez csak idõben ugrik
                        w(j)=w(j)*beta/betaopt; % uj béta korrekció
                        dt=-1/lambdaopt*log(rand); % életidõ sorsolás torzítottan
%                         if dt > tMax && t == 0
%                             s = s + 1;
%                             Egyszer = Egyszer + SumP(k,j);
%                         end
                        t = t + dt;
                        if t > interval(k) % ha az idõ túllépi az adott intervallumot, akkor a következõbe lép
                            w1 = exp(-(interval(k) - (t - dt))*(lambda - lambdaopt)); %(t - dt) helyett t volt
                            w(j) = w(j)*w1; %súlykorrekció
                            t = interval(k); % szakaszhatárra állítás
                            wT(k,j) = w(j); % az adott intervallum határON a neutronlánc súlya
                           
                            if t < max(interval)  %!!!!!! <= nem kell csak <, mert akkor az utolsónál is belép, de ott nem kell
                                %az utsó intervallumnál nem kell
                                dt = -1/lambdaopt*log(rand);
                                t = t + dt; % ha dt2 > resolution (5), k kettõt lép!!
                                w1 = lambda/lambdaopt*exp(-dt*(lambda - lambdaopt));
                                w(j) = w(j)*w1;
                              
                            end
                            if k + 1 + floor(dt/resolution) <= length(interval)%ha nem ugrik ki
                                k = k + 1 + floor(dt/resolution); % egynél többet lép itt esetenként! nem 5tel hanem a lépéssel kell leosztani
                            else % (10,20),30,40 re is ugorhat, ilyenkor 20ra kell állítani, különben k index túllép
                                k = length(interval);
                            end
                        elseif  t < interval(k) % ha nem lépi túl az intervallumhatárt, akkor is kell korrekció
                            
                            w0 = lambda/lambdaopt*exp(-dt*(lambda - lambdaopt));
                            w(j) = w(j)*w0;       
                        end
                    else
                        lefut = lefut + 1;
                        w(j) = w(j)*(1 - beta)/(1 - betaopt); % béta korr.
                    end
                    
                else  % szóródik, súlykorrekció itt is kell
                    w(j) = w(j)*(pS/(pF + pS))/(1 - qopt);
                end
            end
        end
% EgyszerNorm = Egyszer/sum(sum(SumP));
% sumPower = SumPossz;
% sum1JumpPower = Egyszer;
szorzo(z) = lambda*exp(-lambda*(init(z)));
fiplus(z) = sum(sum(SumP)); % össz leadott teljesítmény
end
pdf = szorzo.*fiplus/(szorzo(1)*fiplus(1));
%pdf = szorzo.*fiplus/fiplus(1);
for i = 1:length(interval) % numerikus integrálás trapéz módszerrel
    if i == 1
        integral(i) = (pdf(1) + pdf(1 + i))/2*resolution;
    else
        integral(i) = integral(i - 1) + (pdf(i) + pdf(1 + i))/2*resolution;
    end
end
integral = integral/integral(end);
mcX = zeros(length(integral),1);
nSample = 10000;
for n = 1:nSample
r = rand;
for i = 1:length(integral)
    if i == 1
        if r < integral(1)
            mcX(1) = mcX(1) + 1;
        end
    else
        if r < integral(i) && r >= integral(i - 1)
            mcX(i) = mcX(i) + 1;
        end
    end
end
end
mcX = mcX/mcX(1);

figure(1)
hold on
plot(init,fiplus)
xlabel('Time','fontsize',18)
grid on
%intervallumok határán lévõ súlyok szórása
ylabel('Future power output','fontsize',18)
title({'Time dependence of the future','power output'},'fontsize',16)
hold off

figure(2)
hold on
plot(init,pdf,'b')
plot(init(1:end-1),integral,'r')
plot(init(1:end-1),mcX,'g')
xlabel('Time','fontsize',18)
grid on
%intervallumok határán lévõ súlyok szórása
ylabel({'Weighted future output, its','numerical integral and MC approximation'},'fontsize',14)
legend('Weighted future output (PDF)','Integral (CDF)','MC approximation')
title({'Time dependence of the',' future power output and its integral','and its MC approximation'},'fontsize',14)
hold off



