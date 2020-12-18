function [w,W,wT,SumP,varq,varq2,varq3,lambda,beta,ujlambda,qopt,betaopt,lambdaopt,quotient] = betaOpt(pF,pS,nNeutron,tMax,resolution)
% param�tertorz�t�sos szimul�ci� sz�r�sminimaliz�l�s �letid�sorsol�ssal
%tMax - maximum id�
%resolution - max. id�t h�nyas�val l�ptesse (tMax oszthat� legyen vele)
%iLim = 0 - sima anal�g
%iLim = 1 - szakaszhat�rra vissza�ll��s t�ll�p�s eset�n
warning('off','all')
close all
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
iMax = 40; % q(i), �s q2(i) param�ter f�gg�s�nek felbont�sa
bMax = 40; % �j b�ta param�ter f�gg�s�nek felbont�sa
%ujbeta = 0.3;
if rem(tMax,resolution) == 0
    interval = resolution:resolution:tMax; % id�intervallumok
else %hiba�zenet
    error('A max id� legyen oszthat� a felbont�ssal!');
    exit(0);
end
varq = zeros(bMax,iMax); % relat�v sz�r�sok m�trixa
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
        %minden q(i)-ra m�s:
        wT=zeros(length(interval),nNeutron);
        w = ones(1,nNeutron); % neutronok s�lyai
        W = zeros(length(interval),nNeutron); % s�lyok bankja
        ujbeta(b) =  0.0032 + 0.0001*b;
        for j = 1:nNeutron %neutronl�ptet�s
            t = 0;
            %t = -1/lambda*log(rand); % kezdeti id�t sorsolunk minden neutronnak
            k = 1; % intervallum l�ptet�
            while t < max(interval)
                w(j) = (1 - pA)*w(j); %t�l�l�s vals�ge, implicit capture
                if rand  < q(i) % hasad�s,prompt neutronok csak
                    w(j) = pF/(pF + pS)/q(i)*w(j);
                    %most SumP = fiplus
                    SumP(i,b) = SumP(i,b) + w(j); % leadott �sszteljes�tm�ny (hasad�sok sz�ma, s�lyok �sszege)(i,k) volt, most fiplus
                    W(k,j) = W(k,j) + w(j); % leadott teljes�tm�ny neutrononk�nt, intervallumonk�nt elt�rolva: sum(W(k,:)) = SumP(k) teljes�l
                    w(j) = w(j)*nu;
                    
                    if rand < ujbeta(b) % itt keletkezik k�ss�neutron: ez csak id�ben ugrik
                        w(j)=w(j)*beta/ujbeta(b); % uj b�ta korrekci�
                        dt=-1/ujlambda(i)*log(rand); % �letid� sorsol�s torz�tottan
                        t = t + dt;
                        if t > interval(k) % ha az id� t�ll�pi az adott intervallumot, akkor a k�vetkez�be l�p
                            w1 = exp(-(interval(k) - (t - dt))*(lambda - ujlambda(i))); %(t - dt) helyett t volt
                            w(j) = w(j)*w1; %s�lykorrekci�
                            t = interval(k); % szakaszhat�rra �ll�t�s
                            wT(k,j) = w(j); % az adott intervallum hat�rON a neutronl�nc s�lya
                            maxfel= maxfel + 1;
                            if t < max(interval)  %!!!!!! <= nem kell csak <, mert akkor az utols�n�l is bel�p, de ott nem kell
                                %az uts� intervallumn�l nem kell
                                dt = -1/ujlambda(i)*log(rand);
                                t = t + dt; % ha dt2 > resolution (5), k kett�t l�p!!
                                w1 = lambda/ujlambda(i)*exp(-dt*(lambda - ujlambda(i)));
                                w(j) = w(j)*w1;
                                lefut = lefut + 1;
                            end
                            if k + 1 + floor(dt/resolution) <= length(interval)%ha nem ugrik ki
                                k = k + 1 + floor(dt/resolution); % egyn�l t�bbet l�p itt esetenk�nt! nem 5tel hanem a l�p�ssel kell leosztani
                            else % (10,20),30,40 re is ugorhat, ilyenkor 20ra kell �ll�tani, k�l�nben k index t�ll�p
                                k = length(interval);
                            end
                        elseif  t < interval(k) % ha nem l�pi t�l az intervallumhat�rt, akkor is kell korrekci�
                            
                            w0 = lambda/ujlambda(i)*exp(-dt*(lambda - ujlambda(i)));
                            w(j) = w(j)*w0;
                            bel = bel + 1;
                        end
                    else
                        w(j) = w(j)*(1 - beta)/(1 - ujbeta(b)); % b�ta korr.
                    end
                    
                else  % sz�r�dik, s�lykorrekci� itt is kell
                    w(j) = w(j)*(pS/(pF + pS))/(1 - q(i));
                end
            end
        end
        %for k = 1:length(interval) %sz�r�sok, most csak a v�gs� intervallumot fogja n�zni b�ta' �s lambda' f�ggv�ny�ben
        Ws=W(length(interval),W(length(interval),:)>0); %hossz=length(Ws); %1 volt mindig az 1. k helyett:(
        varq(b,i) = sqrt(sum(Ws.^2)/(sum(Ws)^2) - 1/nNeutron); % 0-k nincsenek benne, r�videbb
        varq2(b,i) = var(wT(length(interval),wT(length(interval),:)>0)); % sima sz�r�s a v�g�n
        Ws3=wT(length(interval),wT(length(interval),:)>0); %hossz=length(Ws3); %relat�v sz�r�s a v�g�n
        varq3(b,i) = sqrt(sum(Ws3.^2)/(sum(Ws3)^2) - 1/nNeutron);
        % end
        
    end
end

figure(1)
hold on
surf(ujlambda,ujbeta,varq)
xlabel('\lambda''','fontsize',18)
grid on
% ahol nincs ott NaN a sz�r�s, a  sz�ls�s�ges param�tertorz�t�s miatt, azaz kihalt
%pontokat kell figyelni, interval nem a 0-t�l indul, az els� pont az 1.
%intervallumban a sz�r�s
ylabel('\beta''','fontsize',18)
zlabel('RSD','fontsize',18)
title({'Biasing parameter dependence of the RSD','of power output in final interval'},'fontsize',16)
hold off

figure(2)
hold on
surf(ujlambda,ujbeta,varq3)
xlabel('\lambda''','fontsize',18)
grid on
%intervallumok hat�r�n l�v� s�lyok sz�r�sa
ylabel('\beta''','fontsize',18)
zlabel('RSD','fontsize',18)
title({'RSD of weight distribution','at the end of the final interval'},'fontsize',16)
hold off

figure(3)
hold on
%surf(ujlambda,ujbeta,SumP') %fiplus nem kell mert sumP m�r az eg�szet
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
%intervallumok hat�r�n l�v� s�lyok sz�r�sa
ylabel('\beta''','fontsize',18)
zlabel('Variance','fontsize',18)
title({'Variance of weight distribution','at the end of the final interval'},'fontsize',16)
hold off
end
