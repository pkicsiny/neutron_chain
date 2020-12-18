function [w,W,wT,SumP,varq,lefut,fiplus] = phibeta( sigmaA,sigmaF,sigmaS,nNeutron,tMax,resolution)
% param�tertorz�t�sos szimul�ci� sz�r�sminimaliz�l�s �letid�sorsol�ssal
%tMax - maximum id�
%resolution - max. id�t h�nyas�val l�ptesse (tMax oszthat� legyen vele)
%iLim = 0 - sima anal�g
%iLim = 1 - szakaszhat�rra vissza�ll��s t�ll�p�s eset�n
warning('off','all')
close all

if ~(tMax > 0 && resolution > 0 && sigmaA >= 0 && sigmaF >= 0 && sigmaS >= 0 && nNeutron > 0)
    error('Nem megfelel� input param�terek!');
    exit(0);
end
betai=[21 142 128 257 75 27]*1e-5;
T12i=[55.7 22.7 6.2 2.3 0.615 0.23];
beta=sum(betai);
lambda=1/( sum(T12i.*betai)/beta );

nu = 2.5; % hasad�si sokszoroz�d�s
iMax = 200; % q(i), �s q2(i) param�ter f�gg�s�nek felbont�sa
sigmaT = sigmaA + sigmaS + sigmaF;
pA=sigmaA/sigmaT;  %vals�gek [0,1]
pF=sigmaF/sigmaT;
pS=sigmaS/sigmaT;
if rem(tMax,resolution) == 0
    interval = resolution:resolution:tMax; % id�intervallumok
else
    error('A max id� legyen oszthat� a felbont�ssal!');
    exit(0);
end

fiplus = zeros(length(interval) + 1,iMax); % j�v�beli teljes�tm�nyek adott t-t�l kezdve: els� sor a teljes, utols� sor 0
fiplus(length(interval) + 1,:) = 0; %utols� 0 pl. 20 a max id� �s 20-t�l ind�tom (csak a sz�ps�g kedv��rt)
init = [0,interval];%(1:end-1)]; [0,5,10,15,20,..]

for z = 1:length(interval) % j�v�beli teljes�tm�ny, sz�zen ind�tom mindig
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
    q(i) = (sigmaT - sigmaS)/sigmaT;
    q2(i) = 0.083625 + 0.00083625*i;
    for j = 1:nNeutron %neutronl�ptet�s
        t = init(z); %% most mindig k�s�bbr�l indul  |__|__|__|__|
        k = z; % intervallum l�ptet� %%1, most mindig k�s�bbr�l indul
        while t < max(interval)
            w(j) = (1 - pA)*w(j); %t�l�l�s vals�ge, implicit capture
            if rand  < q(i) % hasad�s,prompt neutronok csak
                w(j) = pF/(pF + pS)/q(i)*w(j);
                SumP(i,k) = SumP(i,k) + w(j); % leadott �sszteljes�tm�ny (hasad�sok sz�ma, s�lyok �sszege) az adott intervallumban,
                % minden neutron �ltal, mert a t�bbi neutron� is ehhez fog
                % hozz�ad�dni (neutronpopul�ci�)
                W(k,j) = W(k,j) + w(j); % leadott teljes�tm�ny neutrononk�nt, intervallumonk�nt elt�rolva: sum(W(k,:)) = SumP(k) teljes�l
                w(j) = w(j)*nu;
                
                    dt=-1/q2(i)*log(rand); % �letid� sorsol�s torz�tottan
                    t = t + dt;
                    if t > interval(k) % ha az id� t�ll�pi az adott intervallumot, akkor a k�vetkez�be l�p
                        w1 = exp(-(interval(k) - (t - dt))*(lambda - q2(i))); %(t - dt) helyett t volt
                        w(j) = w(j)*w1; %s�lykorrekci�
                        t = interval(k); % szakaszhat�rra �ll�t�s
                        wT(k,j) = w(j); % az adott intervallum hat�rON a neutronl�nc s�lya
                        maxfel= maxfel + 1;
                        if t < max(interval)  %!!!!!! <= nem kell csak <, mert akkor az utols�n�l is bel�p, de ott nem kell
                            %az uts� intervallumn�l nem kell
                            dt = -1/q2(i)*log(rand);
                            t = t + dt; % ha dt2 > resolution (5), k kett�t l�p!!
                            w1 = lambda/q2(i)*exp(-dt*(lambda - q2(i)));
                            w(j) = w(j)*w1;
                            lefut = lefut + 1;
                        end
                        if k + 1 + floor(dt/resolution) <= length(interval)%ha nem ugrik ki
                            k = k + 1 + floor(dt/resolution); % egyn�l t�bbet l�p itt esetenk�nt! nem 5tel hanem a l�p�ssel kell leosztani
                        else % (10,20),30,40 re is ugorhat, ilyenkor 20ra kell �ll�tani, k�l�nben k index t�ll�p
                            k = length(interval);
                        end
                    elseif  t < interval(k) % ha nem l�pi t�l az intervallumhat�rt, akkor is kell korrekci�
                        
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
        Ws=W(k,W(k,:)>0); %hossz=length(Ws); %1 volt mindig az 1. k helyett:(
        varq(k,i) = sum(Ws.^2)/(sum(Ws)^2) - 1/nNeutron; % 0-k nincsenek benne, r�videbb
        varq2(k,i) = var(wT(k,wT(k,:)>0)); % sima sz�r�s a hat�rokon
        Ws3=wT(k,wT(k,:)>0); %hossz=length(Ws3); %hat�ron s�ly rel sz�r�s
        varq3(k,i) = sum(Ws3.^2)/(sum(Ws3)^2) - 1/nNeutron;
    end
end
fiplus(z,:) = sum(SumP'); % �ssz leadott teljes�tm�ny
end

figure(1)
hold on
surf(q2,init,fiplus)
xlabel('\lambda''','fontsize',18)
grid on
%intervallumok hat�r�n l�v� s�lyok sz�r�sa
ylabel('Time interval','fontsize',18)
zlabel('Total power output','fontsize',18)
title({'Biasing parameter and time dependence of the total','power output'},'fontsize',16)
hold off
end
