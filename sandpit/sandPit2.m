function [ sample,integral,monteCarloIntegral,relDev,time,check] = sandPit2( lambda,a,nSample,method,w,i,check)
%g = lambda*exp(-lambda*x);
%f = x/fnorm;
%p = 1/a
tic
if (method == 1 || method == 0) && (w == 1 || w == 0)
if w == 0
    fnorm = a^2/2;
elseif w == 1
    fnorm = 1;
end
integral = (1 - exp(-lambda*a)*(1 + lambda*a))/lambda;%ugyanaz (f*g)
subplot(3,5,i)
hold on
fplot(@(t) t/fnorm,[0,a],'b'); %f plot eredeti
fplot(@(t) t*lambda*exp(-lambda*t),[0,a],'r'); % f*g plot, ami alatti terület kell
if w == 1
fplot(1/a,[0,a],'y');
end
grid on

%módszerek////////////////////////////////////////////////////////////////
if method == 0  %normálásos
    [sample] = calculateSampleNorm(a,nSample,fnorm,w);
elseif method == 1  %rejekciós
    [sample] =calculateSampleRejection(a,nSample,fnorm,w);
end

%relatív szórás és MC integrál////////////////////////////////////////////
if w == 1
    weight = f(sample,fnorm)./p(sample,a);
elseif w == 0
    weight = 1;
end
monteCarloIntegral =fnorm*sum(weight.*g(sample,lambda))/nSample;%mert I/fnorm-ot kapom meg amúgy
relDev = sqrt(sum(((weight.*g(sample,lambda)).^2))/sum(weight.*g(sample,lambda))^2 - 1/nSample);
%plot/////////////////////////////////////////////////////////////////////
[x, fmc] = calculatePDF2(sample, 100, 1); %gyakoriság
%fmc = fmc.*fnorm; % integrál elõtti fnorm, rossz volt az elm. plot,
%f/fnorm kell s akk ez nem kell
subplot(3,5,i)
plot(x,fmc,'g');
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
    fx(i) = sample(i)/fnorm;
end
end

%detektor függvény (g(x))/////////////////////////////////////////////////
function [ gx ] = g( sample,lambda)
for i = 1:length(sample)
    gx(i) = lambda*exp(-lambda*sample(i));%*fnorm kiemelhetõ az egész elé
end
end

%súly/////////////////////////////////////////////////////////////////////
function [px] = p(sample,a)
for i = 1:length(sample)
    px(i) = 1/a; %egyenletes eloszlás pdf
end
end

%rejekciós eljárás////////////////////////////////////////////////////////
function [sample] =calculateSampleRejection(a,nSample,fnorm,w)
sample = zeros(1,nSample);
i = 1;
if w == 0 % nincs súly
    while i<= nSample
        preSample = sqrt(2*rand*fnorm);
        if preSample <= a
            sample(i) = preSample;
            i = i + 1;
        end
    end
elseif w == 1 % van súly
    while i<= nSample
        preSample = c*rand; % Pbõl minták
        if preSample <= a
            sample(i) = preSample;
            i = i + 1;
        end
    end
end
end

%normálásos eljárás///////////////////////////////////////////////////////
function [sample] = calculateSampleNorm(a,nSample,fnorm,w)
if w == 0 % nincs súly
    for i = 1:nSample
        sample(i) = sqrt(rand*2*fnorm); 
    end
elseif w == 1 % van súly
    for i = 1:nSample
        sample(i) = a*rand*fnorm;
    end
end
end
%mintavételezés P bõl ha van súly
