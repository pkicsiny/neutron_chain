function [ sample,integral,monteCarloIntegral,relDev,time,check] = nullDev(a,nSample,method,w,i,check)
%f = x^2/fnorm
%g = x;
%p = (f*g)/integral
tic
if (method == 1 || method == 0) && (w == 1 || w == 0)
c = a; % s�ly pdf egyenletes eloszl�s
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
fplot(@(t) t^3,[0,a],'r'); % f*g plot, ami alatti ter�let kell
if w == 1
fplot(@(t) t^3/integral,[0,a],'y'); % p plot
end
grid on

%m�dszerek////////////////////////////////////////////////////////////////
if method == 0  %norm�l�sos
    [sample] = calculateSampleNorm(nSample,fnorm,w,c);
elseif method == 1  %rejekci�s
    [sample] =calculateSampleRejection(a,nSample,fnorm,w,c);
end

%relat�v sz�r�s �s MC integr�l////////////////////////////////////////////
if w == 1
    weight = f(sample,fnorm)./p(sample,c);
elseif w == 0
    weight = 1;
end
monteCarloIntegral = fnorm*sum(weight.*g(sample))/nSample; %integr�l el�tti fnorm szorz�
relDev = sqrt(sum(((weight.*g(sample)).^2))/sum(weight.*g(sample))^2 - 1/nSample);
%plot/////////////////////////////////////////////////////////////////////
[x, fmc] = calculatePDF2(sample, 100, 1); %gyakoris�g
%fmc = fmc.*fnorm;
subplot(3,5,i)
plot(x,fmc,'g'); %0-kat is hozz�adja
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
ylabel('f(x) gyakoris�g')
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf Mintav�telezett �s elm�leti f(x)','HorizontalAlignment','center','VerticalAlignment', 'top')
grid on
hold off
else
    fprintf('Norm�l�sos: 0, Rejekci�s: 1\nNincs s�ly: 0, Van s�ly: 1\n');
end
time = toc;
end
%f(x) PDF-et ide kell be�rni, hogy mi/////////////////////////////////////
function [fx] = f(sample,fnorm)
for i = 1:length(sample)
    fx(i) = sample(i)^2/fnorm;
end
end

%detektor f�ggv�ny (g(x))/////////////////////////////////////////////////
function [ gx ] = g( sample)
for i = 1:length(sample)
    gx(i) = sample(i); % ide kell �rni g-t
end
end

%s�ly/////////////////////////////////////////////////////////////////////
function [px] = p(sample,c)
for i = 1:length(sample)
    px(i) = sample(i)^3/(c^4/4);
end
end

%rejekci�s elj�r�s////////////////////////////////////////////////////////
function [sample] =calculateSampleRejection(a,nSample,fnorm,w,c)
sample = zeros(1,nSample);
i = 1;
if w == 0 % nincs s�ly
    while i<= nSample
        preSample = nthroot(3*fnorm*rand,3);
        if preSample <= a
            sample(i) = preSample;
            i = i + 1;
        end
    end
elseif w == 1 % van s�ly
    while i<= nSample
        preSample = nthroot(4*c^4/4*rand,4); % Pb�l mint�k
        if preSample <= a
            sample(i) = preSample;
            i = i + 1;
        end
    end
end
end

%norm�l�sos elj�r�s///////////////////////////////////////////////////////
function [sample] = calculateSampleNorm(nSample,fnorm,w,c)
if w == 0 % nincs s�ly
    for i = 1:nSample
        sample(i) = nthroot(3*fnorm*rand,3);
    end
elseif w == 1 % van s�ly
    for i = 1:nSample
        sample(i) = nthroot(4*c^4/4*rand,4);% fnorm=1, ez a p m�r pdf, le van norm�lva
    end
end
end
%mintav�telez�s P b�l ha van s�ly
