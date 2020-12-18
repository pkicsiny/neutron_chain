function [ sample,integral,monteCarloIntegral,relDev,time,check] = nullDev(a,nSample,method,w,i,check)
%f = x^2/fnorm
%g = x;
%p = (f*g)/integral
tic
if (method == 1 || method == 0) && (w == 1 || w == 0)
c = a; % súly pdf egyenletes eloszlás
if w == 0
    fnorm = a^3/3;
elseif w == 1
    fnorm = 1;
end
integral = a^4/4;
%figure(1)
subplot(3,5,i)
hold on
fplot(@(t) t^2/fnorm,[0,a],'b'); %f plot eredeti, max nthroot(3,3)-ig
fplot(@(t) t^3,[0,a],'r'); % f*g plot, ami alatti terület kell
if w == 1
fplot(@(t) t^3/integral,[0,a],'y'); % p plot
end
grid on

%módszerek////////////////////////////////////////////////////////////////
if method == 0  %normálásos
    [sample] = calculateSampleNorm(nSample,fnorm,w,c);
elseif method == 1  %rejekciós
    [sample] =calculateSampleRejection(a,nSample,fnorm,w,c);
end

%relatív szórás és MC integrál////////////////////////////////////////////
if w == 1
    weight = f(sample,fnorm)./p(sample,c);
elseif w == 0
    weight = 1;
end
monteCarloIntegral = fnorm*sum(weight.*g(sample))/nSample; %integrál elõtti fnorm szorzó
relDev = sqrt(sum(((weight.*g(sample)).^2))/sum(weight.*g(sample))^2 - 1/nSample);
%plot/////////////////////////////////////////////////////////////////////
[x, fmc] = calculatePDF2(sample, 100, 1); %gyakoriság
%fmc = fmc.*fnorm;
subplot(3,5,i)
plot(x,fmc,'g'); %0-kat is hozzáadja
str = sprintf('%d db minta', nSample);
title(str);
if check == 0
if w == 1
legend('f(x)','f(x)*g(x)','p(x)','MC PDF')
check = 1;
else
legend('f(x)','f(x)*g(x)','MC PDF')
check = 1;
end
end
xlabel('x')
ylabel('f(x) gyakoriság')
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf Mintavételezett és elméleti f(x)','HorizontalAlignment','center','VerticalAlignment', 'top')
grid on
hold off
else
    fprintf('Normálásos: 0, Rejekciós: 1\nNincs súly: 0, Van súly: 1\n');
end
time = toc;
end
%f(x) PDF-et ide kell beírni, hogy mi/////////////////////////////////////
function [fx] = f(sample,fnorm)
for i = 1:length(sample)
    fx(i) = sample(i)^2/fnorm;
end
end

%detektor függvény (g(x))/////////////////////////////////////////////////
function [ gx ] = g( sample)
for i = 1:length(sample)
    gx(i) = sample(i); % ide kell írni g-t
end
end

%súly/////////////////////////////////////////////////////////////////////
function [px] = p(sample,c)
for i = 1:length(sample)
    px(i) = sample(i)^3/(c^4/4);
end
end

%rejekciós eljárás////////////////////////////////////////////////////////
function [sample] =calculateSampleRejection(a,nSample,fnorm,w,c)
sample = zeros(1,nSample);
i = 1;
if w == 0 % nincs súly
    while i<= nSample
        preSample = nthroot(3*fnorm*rand,3);
        if preSample <= a
            sample(i) = preSample;
            i = i + 1;
        end
    end
elseif w == 1 % van súly
    while i<= nSample
        preSample = nthroot(4*c^4/4*rand,4); % Pbõl minták
        if preSample <= a
            sample(i) = preSample;
            i = i + 1;
        end
    end
end
end

%normálásos eljárás///////////////////////////////////////////////////////
function [sample] = calculateSampleNorm(nSample,fnorm,w,c)
if w == 0 % nincs súly
    for i = 1:nSample
        sample(i) = nthroot(3*fnorm*rand,3);
    end
elseif w == 1 % van súly
    for i = 1:nSample
        sample(i) = nthroot(4*c^4/4*rand,4);% fnorm=1, ez a p már pdf, le van normálva
    end
end
end
%mintavételezés P bõl ha van súly
