warning('off','all')
close all
clear all
pF = 0.25;
pS = 0.25;
nNeutron = 10000;
tMax = 10;
resolution = 5;
iMax = 10;
pA = 1 - pF - pS;
if ~(tMax > 0 && resolution > 0 && pA >= 0 && pF >= 0 && pS >= 0 && nNeutron > 0 && pF <= 1 && pS <= 1)
    error('Nem megfelel� input param�terek!');
    exit(0);
end
%k�rd�s majd??????
betai=[21 142 128 257 75 27]*1e-5;
T12i=[55.7 22.7 6.2 2.3 0.615 0.23];
beta=sum(betai);
lambda=1/( sum(T12i.*betai)/beta );

nu = 2.5; % hasad�si sokszoroz�d�s
if rem(tMax,resolution) == 0
    interval = resolution:resolution:tMax; % id�intervallumok
else %hiba�zenet
    error('A max id� legyen oszthat� a felbont�ssal!');
    exit(0);
end
qopt = 1 - pS;
betaopt = 1 - pF*nu*(1 - beta)/qopt;
lambdaopt = lambda*pF*nu*beta/(qopt*betaopt);
quotientprompt = pS + (1 - beta)*pF*nu;
quotientall = pS + pF*nu;
init = [0,interval];%(1:end-1)]; [0,5,10,15,20,..]
Jump1Chain = zeros(iMax,1);
Jump2Chain = zeros(iMax,1);
%for z = 1:length(interval)
for i = 1:iMax
varq = zeros(length(interval),1); % relat�v sz�r�sok m�trixa
SumP = zeros(length(interval),nNeutron);
SumPossz = 0;
Egyszer = 0;
Ketszer = 0;
lefut = 0;
        wT=zeros(length(interval),nNeutron);
        w = ones(1,nNeutron); % neutronok s�lyai
        W = zeros(length(interval),nNeutron); % s�lyok bankja
        for j = 1:nNeutron %neutronl�ptet�s
            step = 0; %h�nyat ugrik mire kirep�l
            %t = init(z);
            t = 0;
            %k = z; % intervallum l�ptet�
           k = 1;
            while t < max(interval)
              rr=rand;
              w(j) = (1 - pA)*w(j);
                if rr  < pF/(1-pA) % hasad�s,prompt neutronok csak
                   SumP(k,j) = SumP(k,j) + w(j);
                    SumPossz = SumPossz + w(j); % leadott �sszteljes�tm�ny (hasad�sok sz�ma, s�lyok �sszege)(i,k) volt
                    W(k,j) = W(k,j) + w(j); % leadott teljes�tm�ny neutrononk�nt, intervallumonk�nt elt�rolva: sum(W(k,:)) = SumP(k) teljes�l
                    w(j) = w(j)*nu;
                    r = rand;
                    if r < beta % itt keletkezik k�ss�neutron: ez csak id�ben ugrik
                        dt=-1/lambda*log(rand); % �letid� sorsol�s torz�tottan
                        step = step + 1;
                        if dt > tMax && step == 1
                            Jump1Chain(i) = Jump1Chain(i) + 1;
                            Egyszer = Egyszer + SumP(k,j);
                        elseif dt > (tMax - t) && step == 2
                            Jump2Chain(i) = Jump2Chain(i) + 1;
                            Ketszer = Ketszer + sum(SumP(:,j));
                        end
                        t = t + dt;
                        if t > interval(k) % ha az id� t�ll�pi az adott intervallumot, akkor a k�vetkez�be l�p
                            t = interval(k); % szakaszhat�rra �ll�t�s
                            wT(k,j) = w(j); % az adott intervallum hat�rON a neutronl�nc s�lya
                            if t < max(interval)  %!!!!!! <= nem kell csak <, mert akkor az utols�n�l is bel�p, de ott nem kell
                                dt = -1/lambda*log(rand);
                                t = t + dt;   
                            end
                            if k + 1 + floor(dt/resolution) <= length(interval)%ha nem ugrik ki
                                k = k + 1 + floor(dt/resolution); % egyn�l t�bbet l�p itt esetenk�nt! nem 5tel hanem a l�p�ssel kell leosztani
                            else % (10,20),30,40 re is ugorhat, ilyenkor 20ra kell �ll�tani, k�l�nben k index t�ll�p
                                k = length(interval);
                            end
                        elseif  t < interval(k) % ha nem l�pi t�l az intervallumhat�rt, akkor is kell korrekci�
                            
                            %w0 = lambda/lambdaopt*exp(-dt*(lambda - lambdaopt));
                            %w(j) = w(j)*w0;
                          
                        end
                    else
                        lefut = lefut + 1;
                   
                    end
                end
            end
        end
EgyszerNorm(i) = Egyszer/SumPossz;
KetszerNorm(i) = Ketszer/SumPossz;
sumPower(i) = SumPossz;
sum1JumpPower(i) = Egyszer;
sum2JumpPower(i) = Ketszer;
end
Avg1JumpPowerNormed = sum(EgyszerNorm)/iMax;
Avg2JumpPowerNormed = sum(KetszerNorm)/iMax;
Avg1JumpChainNormed = sum(Jump1Chain)/iMax/nNeutron;
Avg2JumpChainNormed = sum(Jump2Chain)/iMax/nNeutron;
AvgSumP = sum(sumPower)/iMax;
Avg1JumpPower = sum(sum1JumpPower)/iMax;
Avg2JumpPower = sum(sum2JumpPower)/iMax;
analyticalQuotient = pF*nu + pS;
analyticalSumP = nNeutron*pF/(1 - (pS + nu*pF*(1 - beta)));


