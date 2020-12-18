function [w,W,SumP,varq,lefut] = delayedNeutronChain( sigmaA,sigmaF,sigmaS,nNeutron,tMax,resolution,iLim)
    % param�tertorz�t�sos szimul�ci� sz�r�sminimaliz�l�s �letid�sorsol�ssal
    %tMax - maximum id�
    %resolution - max. id�t h�nyas�val l�ptesse (tMax oszthat� legyen vele)
    %iLim = 0 - sima anal�g
    %iLim = 1 - szakaszhat�rra vissza�ll��s t�ll�p�s eset�n
warning('off','all')
close all

if ~(tMax > 0 && resolution > 0 && sigmaA >= 0 && sigmaF >= 0 && sigmaS >= 0 && nNeutron > 0 && (iLim == 0 || iLim == 1 || iLim == 2))
error('Nem megfelel� input param�terek!');
exit(0);
end
%k�rd�s majd??????
betai=[21 142 128 257 75 27]*1e-5;
T12i=[55.7 22.7 6.2 2.3 0.615 0.23];
beta=sum(betai);
lambda=1/( sum(T12i.*betai)/beta );

nu = 2.5; % hasad�si sokszoroz�d�s
iMax = 300; % q(i), �s q2(i) param�ter f�gg�s�nek felbont�sa
sigmaT = sigmaA + sigmaS + sigmaF;
pA=sigmaA/sigmaT;  %vals�gek [0,1]
pF=sigmaF/sigmaT;
pS=sigmaS/sigmaT;
if rem(tMax,resolution) == 0
interval = [resolution:resolution:tMax]; % id�intervallumok
else
    error('A max id� legyen oszthat� a felbont�ssal!');
    exit(0);
end
varq = zeros(length(interval),iMax); % relat�v sz�r�sok m�trixa
SumP = zeros(iMax,length(interval));
lefut = 0;
bel = 0;
maxfel = 0;
for i = 1:iMax
    %minden q(i)-ra m�s:
    wT=zeros(length(interval),nNeutron);
    w = ones(1,nNeutron); % neutronok s�lyai
    W = zeros(length(interval),nNeutron); % s�lyok bankja
    if iLim == 2
        q(i) = (sigmaT - sigmaS)/sigmaT;
        q2(i) = 0.1 + 0.0016*i;
    else
        q(i) = 0.1 + 0.002*i;
    end
    for j = 1:nNeutron %neutronl�ptet�s
        t = 0;
        k = 1; % intervallum l�ptet� 
        while t < max(interval)
            w(j) = (1 - pA)*w(j); %t�l�l�s vals�ge, implicit capture
            if rand  < q(i) % hasad�s,prompt neutronok csak
                w(j) = pF/(pF + pS)/q(i)*w(j); 
                SumP(i,k) = SumP(i,k) + w(j); % leadott �sszteljes�tm�ny (hasad�sok sz�ma, s�lyok �sszege) az adott intervallumban,
                % minden neutron �ltal, mert a t�bbi neutron� is ehhez fog
                % hozz�ad�dni (neutronpopul�ci�)
                W(k,j) = W(k,j) + w(j); % leadott teljes�tm�ny neutrononk�nt, intervallumonk�nt elt�rolva: sum(W(k,:)) = SumP(k) teljes�l
                w(j) = w(j)*nu;
                dt=-1/lambda*log(rand); % �letid� sorsol�s (sigmaT vagy lambda???????)
                if iLim == 2
                dt=-1/q2(i)*log(rand); % �letid� sorsol�s torz�tottan
                end
                t = t + dt;
                if t > interval(k) %&& ter < interval(k)  % ha az id� t�ll�pi az adott intervallumot, akkor a k�vetkez�be l�p
                    if iLim == 1 % szakaszhat�rra vissza�ll�tja az id�t
                        t = interval(k);
                        dt = -1/lambda*log(rand); % �j �letid� sorsol�s dt2
                        t = t + dt;
                    elseif iLim == 2 % vissza�ll�t + visszatorz�t
                        w1 = exp(-(interval(k) - (t - dt))*(lambda - q2(i))); %(t - dt) helyett t volt
                        w(j) = w(j)*w1; %s�lykorrekci�
                        t = interval(k); % szakaszhat�rra �ll�t�s
                        wT(k,j) = w(j); % elmentj�k
                        maxfel= maxfel + 1;
                        if t < max(interval)  %!!!!!! <= nem kell csak <, mert akkor az utols�n�l is bel�p, de ott nem kell
                         %az uts� intervallumn�l nem kell   
                        dt = -1/q2(i)*log(rand);
                        t = t + dt; % ha dt2 > resolution (5), k kett�t l�p!!
                        w1 = lambda/q2(i)*exp(-dt*(lambda - q2(i)));
                        w(j) = w(j)*w1;
                        lefut = lefut + 1;
                        end
                    end
                    if k + 1 + floor(dt/resolution) <= length(interval)%ha nem ugrik ki
                    k = k + 1 + floor(dt/resolution); % egyn�l t�bbet l�p itt esetenk�nt! nem 5tel hanem a l�p�ssel kell leosztani
                    else % (10,20),30,40 re is ugorhat, ilyenkor 20ra kell �ll�tani, k�l�nben k index t�ll�p
                        k = length(interval);
                    end
                    elseif iLim == 2 && t < interval(k) % ha nem l�pi t�l az intervallumhat�rt, akkor is kell korrekci�
                   
                    w0 = lambda/q2(i)*exp(-dt*(lambda - q2(i)));
                    w(j) = w(j)*w0;
                    bel = bel + 1;
                end
            else  % sz�r�dik, s�lykorrekci� itt is kell
                w(j) = w(j)*(pS/(pF + pS))/(1 - q(i));
            end
        end
    end
    for k = 1:length(interval)
        Ws=W(k,W(k,:)>0); hossz=length(Ws); %1 volt mindig az 1. k helyett:(
        varq(k,i) = sqrt(sum(Ws.^2)/(sum(Ws)^2) - 1/nNeutron); % 0-k nincsenek benne, r�videbb
    % varq(k,i) = sqrt(sum(W(k,:).^2)/(sum(W(k,:))^2) - 1/nNeutron); %s�lyok sz�r�sn�gyzete a k. intervallumban q(i) param�ter eset�n
    end
end

figure(1)
hold on
for k = 1:length(interval)
    if iLim == 2
plot(q2,varq(k,:),'.-') %s�lyok sz�r�sa
xlabel('q2 parameter')
    else
        plot(q,varq(k,:),'.-') %s�lyok sz�r�sa
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
% ahol nincs ott NaN a sz�r�s, a  sz�ls�s�ges param�tertorz�t�s miatt, azaz kihalt
%pontokat kell figyelni, interval nem a 0-t�l indul, az els� pont az 1.
%intervallumban a sz�r�s
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