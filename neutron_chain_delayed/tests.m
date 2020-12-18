function [w,W,SumP,varq,lefut] = csak3FPteszt( sigmaA,sigmaF,sigmaS,nNeutron,tMax,resolution)
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

nu = 2.5; 
iMax = 1; 
sigmaT = sigmaA + sigmaS + sigmaF;
pA=sigmaA/sigmaT;  
pF=sigmaF/sigmaT;
pS=sigmaS/sigmaT;
if rem(tMax,resolution) == 0
interval = [resolution:resolution:tMax]; 
else
    error('A max idõ legyen osztható a felbontással!');
    exit(0);
end
varq = zeros(length(interval),iMax); 
SumP = zeros(iMax,length(interval));
lefut = 0;
bel = 0;
maxfel = 0;
lefut2 = 0;
bel2 = 0;
maxfel2 = 0;

for i = 1:iMax 
    w = ones(1,nNeutron);
    W = zeros(length(interval),nNeutron);
    q(i) = (sigmaT - sigmaS)/sigmaT;
    q2(i) = 0.18;
    %q2(i) = 0.1 + 0.002*i;
    wT=zeros(length(interval),nNeutron);

w2=ones(1,nNeutron);SumP12=0.0;SumP22=0.0;SumP32=0.0;SumP42=0.0;W2=zeros(4,nNeutron);
T1=5.0;
T2=10.0;
T3=15.0;
T4=20.0;
wT2=zeros(4,nNeutron);

    for j = 1:nNeutron 
        t = 0;
        k = 1;
        z = 1;
        s = 1;
        while t < max(interval)
            w(j) = (1 - pA)*w(j); 
%%%%
            if rand  < 1%q(i) 
                w(j) = pF/(pF + pS)/q(i)*w(j); 
                SumP(i,k) = SumP(i,k) + w(j); 
                W(k,j) = W(k,j) + w(j); 
                w(j) = w(j)*nu;
                dt(z)=-1/q2(i)*log(rand); 
                t = t + dt(z);
                if t > interval(k) 
                        w1 = exp(-(interval(k) - (t - dt(z)))*(lambda - q2(i))); 
                        w(j) = w(j)*w1; 
                        t = interval(k); 
                        wT(k,j) = w(j);
                        maxfel= maxfel + 1;
                        if t < max(interval)  
                        dt2(s) = -1/q2(i)*log(rand);
                        t = t + dt2(s); 
                        w1 = lambda/q2(i)*exp(-dt2(s)*(lambda - q2(i)));
                        w(j) = w(j)*w1;
                        
                        z = z + 1;
                        s = s + 1;
                        lefut = lefut + 1;
                        end
                    k = k + 1;
                elseif t < interval(k) 
                    w0 = lambda/q2(i)*exp(-dt(z)*(lambda - q2(i)));
                    w(j) = w(j)*w0;
                    
                    z = z + 1;
                    bel = bel + 1;
                end
%%%%
            else  
                w(j) = w(j)*(pS/(pF + pS))/(1 - q(i));
            end
        end
        puska = 1;
        %%
        if puska == 1
            t=0.0;
            z = 1;
            s = 1;
        while t<T4
            w2(1,j)=w2(1,j)*(1-sigmaA/sigmaT); % implicit befogas
            if rand < 1%q(i) % hasadas
                w2(1,j)=w2(1,j)*(sigmaF/(sigmaT-sigmaA)/q(i));
                if(t<T1)
                    SumP12=SumP12+w2(1,j);
                    W2(1,j)=W2(1,j)+w2(1,j);
                elseif(t <T2)
                    SumP22=SumP22+w2(1,j);
                    W2(2,j)=W2(2,j)+w2(1,j);
                elseif(t <T3)
                    SumP32=SumP32+w2(1,j);
                    W2(3,j)=W2(3,j)+w2(1,j);
                elseif(t <T4)
                    SumP42=SumP42+w2(1,j);
                    W2(4,j)=W2(4,j)+w2(1,j);
                end
                w2(1,j)=w2(1,j)*nu;
               %dt=-1.0/q2(i)*log(1.0-rand); % sorsoljunk prekurzor idot
               %dt(z)

                if(t < T1 && t+dt(z) > T1)
                   w12=exp(-(T1-t)*(lambda-q2(i)));
                   w2(1,j)=w2(1,j)*w12;
                   t=T1;
                   wT2(1,j)=w2(1,j);

                  % dt=-1.0/q2(i)*log(1.0-rand);
                   t=t+dt2(s);
                   w12=lambda/q2(i)*exp(-dt2(s)*(lambda-q2(i)));
                   w2(1,j)=w2(1,j)*w12;
                   z = z + 1;
                   s = s + 1;
                   lefut2 = lefut2 + 1;
                   maxfel2= maxfel2 + 1;
                elseif(t < T2 && t+dt(z) > T2)
                    w12=exp(-(T2-t)*(lambda-q2(i)));
                   w2(1,j)=w2(1,j)*w12;
                   t=T2;
                   wT2(2,j)=w2(1,j);

                  % dt=-1.0/q2(i)*log(1.0-rand);
                   t=t+dt2(s);
                   w12=lambda/q2(i)*exp(-dt2(s)*(lambda-q2(i)));
                   w2(1,j)=w2(1,j)*w12;
                   z = z + 1;
                   s = s + 1;
                   lefut2 = lefut2 + 1;
                   maxfel2 = maxfel2 + 1;
               elseif(t < T3 && t+dt(z) > T3)
                   w12=exp(-(T3-t)*(lambda-q2(i)));
                   w2(1,j)=w2(1,j)*w12;
                   t=T3;
                   wT2(3,j)=w2(1,j);

                   %dt=-1.0/q2(i)*log(1.0-rand);
                   t=t+dt2(s);
                   w12=lambda/q2(i)*exp(-dt2(s)*(lambda-q2(i)));
                   w2(1,j)=w2(1,j)*w12;
                   z = z + 1;
                   s = s + 1;
                   lefut2 = lefut2 + 1;
                   maxfel2 = maxfel2 + 1;
                elseif(t < T4 && t+dt(z) > T4)
                   
                   w12=exp(-(T4-t)*(lambda-q2(i)));
                   w2(1,j)=w2(1,j)*w12;
                   wT2(4,j)=w2(1,j);
                   t=T4;
                   maxfel2 = maxfel2 + 1;
                else
                   
                   w02=lambda/q2(i)*exp(-dt(z)*(lambda-q2(i)));
                   w2(1,j)=w2(1,j)*w02;
                   t=t+dt(z);
                   
                   z = z + 1;
                   bel2 = bel2 + 1;
                end
 
            else % szoras
                w2(1,j)=w2(1,j)*(sigmaS/(sigmaT-sigmaA)/(1-q(i)));
            end
        end
        end
        %%
    end
    
    for b = 1:length(interval)
        Ws=W(b,W(b,:)>0); hossz=length(Ws); 
        varq(b,i) = sqrt(sum(Ws.^2)/(sum(Ws)^2) - 1/nNeutron); 
    end
end

figure(1)
hold on
for k = 1:length(interval) 
plot(q2,varq(k,:),'.-') 
xlabel('q2 parameter')
end
grid on
ylabel('Relative deviation')
title('Relative deviation - q2')


%%%%%%%%%%%%%%%%%
 

figure(2)
hold on
grid on
   
    Ws2=W2(1,W2(1,:)>0); N_1_3=length(Ws2);
    var1_2(s)=sqrt(sum(Ws2.*Ws2)/(sum(Ws2)*sum(Ws2))-1/nNeutron); %N_1_3
     
    Ws2=W2(2,W2(2,:)>0); N_2_3=length(Ws2);
    var2_2(s)=sqrt(sum(Ws2.*Ws2)/(sum(Ws2)*sum(Ws2))-1/nNeutron); %N_2_3
     
    Ws2=W2(3,W2(3,:)>0); N_3_3=length(Ws2);
    var3_2(s)=sqrt(sum(Ws2.*Ws2)/(sum(Ws2)*sum(Ws2))-1/nNeutron); %N_3_3
     
    Ws2=W2(4,W2(4,:)>0); N_4_3=length(Ws2);
    var4_2(s)=sqrt(sum(Ws2.*Ws2)/(sum(Ws2)*sum(Ws2))-1/nNeutron); %N_4_3
     

    plot(q2,var1_2,'r.-');
    plot(q2,var2_2,'g.-');
    plot(q2,var3_2,'b.-');
    plot(q2,var4_2,'k.-');
     
    fprintf('nonAnalog, q2= %f: P1=%f \t P2=%f \t P3= %f \t P3= %f\n',q2,SumP12,SumP22,SumP32,SumP42);
    % vissza az alapertekekre
    w2=ones(1,N);SumP12=0.0;SumP22=0.0;SumP32=0.0;SumP42=0.0;W2=zeros(4,N);wT2=zeros(4,N);
    plot(q2,var1_2,'r.-');
    plot(q2,var2_2,'g.-');
    plot(q2,var3_2,'b.-');
    plot(q2,var4_2,'k.-');
legend('T_0 <= t < T_1','T_1 <= t < T_2','T_2 <= t < T_3','T_3 <= t < T_4');
xlabel('\lambda''');
ylabel('\sigma');
hold off
end