warning('off','all')
close all
pF = 0.298;
pS = 0.256;
nNeutron = 1000;
tMax = 100;
resolution = 1;
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

quotient = pS + (1 - beta)*pF*nu;
szorzo = zeros(length(interval) + 1,1);
integral = zeros(length(interval),1);
fiplus = zeros(length(interval) + 1,1); % j�v�beli teljes�tm�nyek adott t-t�l kezdve: els� sor a teljes, utols� sor 0
fiplus(length(interval) + 1) = 0; %utols� 0 pl. 20 a max id� �s 20-t�l ind�tom (csak a sz�ps�g kedv��rt)
init = [0,interval];%(1:end-1)]; [0,5,10,15,20,..]
%SumPossz = zeros(length(interval),1);
s = 0;
for z = 1:length(interval)
varq = zeros(length(interval),1); % relat�v sz�r�sok m�trixa
SumP = zeros(length(interval),nNeutron);
SumPossz = 0;
Egyszer = 0;
lefut = 0;
        wT=zeros(length(interval),nNeutron);
        w = ones(1,nNeutron); % neutronok s�lyai
        W = zeros(length(interval),nNeutron); % s�lyok bankja
        for j = 1:nNeutron %neutronl�ptet�s
            t = init(z);
           % t = 0;
            s2 = 0;
            k = z; % intervallum l�ptet�
           %k = 1;
            while t < max(interval)
                w(j) = (1 - pA)*w(j); %t�l�l�s vals�ge, implicit capture
                if rand  < qopt % hasad�s,prompt neutronok csak
                    w(j) = pF/(pF + pS)/qopt*w(j);
                    SumP(k,j) = SumP(k,j) + w(j);
                    SumPossz = SumPossz + w(j); % leadott �sszteljes�tm�ny (hasad�sok sz�ma, s�lyok �sszege)(i,k) volt
                    W(k,j) = W(k,j) + w(j); % leadott teljes�tm�ny neutrononk�nt, intervallumonk�nt elt�rolva: sum(W(k,:)) = SumP(k) teljes�l
                    w(j) = w(j)*nu;
                    r = rand;
                    if r < betaopt % itt keletkezik k�ss�neutron: ez csak id�ben ugrik
                        w(j)=w(j)*beta/betaopt; % uj b�ta korrekci�
                        dt=-1/lambdaopt*log(rand); % �letid� sorsol�s torz�tottan
%                         if dt > tMax && t == 0
%                             s = s + 1;
%                             Egyszer = Egyszer + SumP(k,j);
%                         end
                        t = t + dt;
                        if t > interval(k) % ha az id� t�ll�pi az adott intervallumot, akkor a k�vetkez�be l�p
                            w1 = exp(-(interval(k) - (t - dt))*(lambda - lambdaopt)); %(t - dt) helyett t volt
                            w(j) = w(j)*w1; %s�lykorrekci�
                            t = interval(k); % szakaszhat�rra �ll�t�s
                            wT(k,j) = w(j); % az adott intervallum hat�rON a neutronl�nc s�lya
                           
                            if t < max(interval)  %!!!!!! <= nem kell csak <, mert akkor az utols�n�l is bel�p, de ott nem kell
                                %az uts� intervallumn�l nem kell
                                dt = -1/lambdaopt*log(rand);
                                t = t + dt; % ha dt2 > resolution (5), k kett�t l�p!!
                                w1 = lambda/lambdaopt*exp(-dt*(lambda - lambdaopt));
                                w(j) = w(j)*w1;
                              
                            end
                            if k + 1 + floor(dt/resolution) <= length(interval)%ha nem ugrik ki
                                k = k + 1 + floor(dt/resolution); % egyn�l t�bbet l�p itt esetenk�nt! nem 5tel hanem a l�p�ssel kell leosztani
                            else % (10,20),30,40 re is ugorhat, ilyenkor 20ra kell �ll�tani, k�l�nben k index t�ll�p
                                k = length(interval);
                            end
                        elseif  t < interval(k) % ha nem l�pi t�l az intervallumhat�rt, akkor is kell korrekci�
                            
                            w0 = lambda/lambdaopt*exp(-dt*(lambda - lambdaopt));
                            w(j) = w(j)*w0;       
                        end
                    else
                        lefut = lefut + 1;
                        w(j) = w(j)*(1 - beta)/(1 - betaopt); % b�ta korr.
                    end
                    
                else  % sz�r�dik, s�lykorrekci� itt is kell
                    w(j) = w(j)*(pS/(pF + pS))/(1 - qopt);
                end
            end
        end
% EgyszerNorm = Egyszer/sum(sum(SumP));
% sumPower = SumPossz;
% sum1JumpPower = Egyszer;
szorzo(z) = lambda*exp(-lambda*(init(z)));
fiplus(z) = sum(sum(SumP)); % �ssz leadott teljes�tm�ny
end
pdf = szorzo.*fiplus/(szorzo(1)*fiplus(1));
%pdf = szorzo.*fiplus/fiplus(1);
for i = 1:length(interval) % numerikus integr�l�s trap�z m�dszerrel
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
%intervallumok hat�r�n l�v� s�lyok sz�r�sa
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
%intervallumok hat�r�n l�v� s�lyok sz�r�sa
ylabel({'Weighted future output, its','numerical integral and MC approximation'},'fontsize',14)
legend('Weighted future output (PDF)','Integral (CDF)','MC approximation')
title({'Time dependence of the',' future power output and its integral','and its MC approximation'},'fontsize',14)
hold off



