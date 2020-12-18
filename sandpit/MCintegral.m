function [integral,monteCarloIntegral,relDev,time,check] = MCintegral(lambda,a,nSample,method,task,i,check)
tic
if (method == 1 || method == 0) && (task == 4 || task == 3 || task == 2 || task == 1)
    
    integral = (1 - exp(-lambda*a)*(1 + lambda*a))/lambda;
    sqrtfnorm = 1/sqrt(lambda*integral)*(1 - exp(-lambda*a));
    if task == 1    % f = lambda*exp(-lambda*x)/fnorm, g = x
        fnorm = 1 - exp(-lambda*a);
    elseif task == 2    % f = x/fnorm, g = lambda*exp(-lambda*x)
        fnorm = a^2/2;
    elseif task == 3    % f = lambda*exp(-lambda*x), g = x, p= 1/a
        fnorm = 1;
    elseif task == 4    % f = lambda*exp(-lambda*x), g = x, p = (f*g)/integral
        fnorm = integral; % a-ig intergálás

            fnorm2a = (1 - exp(-lambda*2*a)*(1 + lambda*2*a))/lambda;
            %2a-ig kell, mert a mintavtelezés is addig van

    end
    
    if i == 1 %kellõ plotok beállítása
        i2 = 1;
    elseif i == 2
        i2 = 2;
    elseif i == 6
        i2 = 3;
    elseif i == 13
        i2 = 4;
    end
    if i == 1 || i == 2 || i == 6 || i == 13
    figure(task)
    subplot(2,2,i2)
    hold on
    if task == 2
        
        fplot(@(t) t/fnorm,[0,a],'m'); %f = x plot eredeti, ami máskor g de most ezt mintavételezzük ezért ennek kell pdfnek lennie
        fplot(@(t) lambda*exp(-lambda*t),[0,a],'b'); %g = lambda*exp(-lambda*x) plot
    else
        fplot(@(t) lambda*exp(-lambda*t)/fnorm,[0,a],'b'); %f = lambda*exp(-lambda*x)/fnorm plot eredeti
        fplot(@(t) t,[0,a],'m'); %g = x plot
    end
    fplot(@(t) t*lambda*exp(-lambda*t),[0,a],'r'); % f*g plot, ami alatti terület kell
    
    if task == 3 % p plot
        fplot(1/a,[0,a],'y');
    elseif task == 4
        fplot(@(t) t*lambda*exp(-lambda*t)/integral,[0,a],'y');
    end
    grid on
    end
    %módszerek////////////////////////////////////////////////////////////////
    if method == 0  %normálásos
        [sample] = calculateSampleNorm(lambda,a,nSample,fnorm,sqrtfnorm,task);
    elseif method == 1  %rejekciós
        [sample] =calculateSampleRejection(lambda,a,nSample,fnorm,sqrtfnorm,task);
    end
    
    %relatív szórás és MC integrál////////////////////////////////////////////
    if task < 3
        weight = 1;
    else
        weight = f(sample,lambda,fnorm,task)./p(sample,lambda,a,integral,task);
    end
    monteCarloIntegral = fnorm*sum(weight.*g(sample,lambda,task))/nSample; %integrál elõtti fnorm szorzó
    if task == 4
        monteCarloIntegral = sum(weight.*g(sample,lambda,task))/nSample; %integrál elõtti fnorm szorzó
    end
    relDev = real(sqrt(sum(((weight.*g(sample,lambda,task)).^2))/sum(weight.*g(sample,lambda,task))^2 - 1/nSample));
    if relDev < 1e-8
        relDev = 0;
    end
    %plot/////////////////////////////////////////////////////////////////////
    [x, fmc] = calculatePDF2(sample, 100, 1); %gyakoriság
    if task == 4
    %fmc = fmc.*sqrtfnorm;
    end
    if i == 1 || i == 2 || i == 6 || i == 13
    figure(task)
    subplot(2,2,i2)
    plot(x,fmc,'g');
    str = sprintf('%d samples', nSample);
    title(str);
    if check == 0
        if task > 2
            legend('f(x)','g(x)','f(x)*g(x)','p(x)','MC PDF')
            check = 1;
        else
            legend('f(x)','g(x)','f(x)*g(x)','MC PDF')
            check = 1;
        end
    end
    xlabel('x','fontsize',16)
    ylabel('Examined functions','fontsize',16)
    set(legend,'FontSize',12);

    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    if task == 1
        text(0.5, 1,'\bf MC integration plot, f(x)=lambda*exp(-lambda*x), g(x)=x','HorizontalAlignment','center','VerticalAlignment', 'top','fontsize',14)
    elseif task == 2
        text(0.5, 1,'\bf Reversed MC integration: f(x)=x, g(x)=lambda*exp(-lambda*x)','HorizontalAlignment','center','VerticalAlignment', 'top','fontsize',14)
    elseif task == 3
        text(0.5, 1,'\bf MC integration with uniform weight function, f(x)=lambda*exp(-lambda*x), g(x)=x, p(x)=1/a','HorizontalAlignment','center','VerticalAlignment', 'top','fontsize',14)
    elseif task == 4
        text(0.5, 1,'\bf Zero standard deviation MC integration, f(x)=lambda*exp(-lambda*x), g(x)=x, p(x)=(f(x)*g(x))/I','HorizontalAlignment','center','VerticalAlignment', 'top','fontsize',14)
    end
    grid on
    hold off
    end
else
    fprintf('Normálásos: 0, Rejekciós: 1\nNincs súly: 0, Van súly: 1\n');
end
time = toc;
end
%f(x) PDF-et ide kell beírni, hogy mi/////////////////////////////////////
function [fx] = f(sample,lambda,fnorm,task)
for i = 1:length(sample)
    if task == 2
        fx(i) = sample(i)/fnorm;
    else
        fx(i) = lambda*exp(-lambda*sample(i));
    end
end
end

%detektor függvény (g(x))/////////////////////////////////////////////////
function [ gx ] = g(sample,lambda,task)
for i = 1:length(sample)
    if task == 2
        gx(i) = lambda*exp(-lambda*sample(i));
    else
        gx(i) = sample(i); % ide kell írni g-t
    end
end
end

%súly/////////////////////////////////////////////////////////////////////
function [px] = p(sample,lambda,a,integral,task)
for i = 1:length(sample)
    if task == 3
        px(i) = 1/a; %egyenletes eloszlás pdf
    elseif task == 4
        px(i) = lambda*exp(-lambda*sample(i))*sample(i)/integral;
    end
end
end

%rejekciós eljárás (1)////////////////////////////////////////////////////////
function [sample] =calculateSampleRejection(lambda,a,nSample,fnorm,sqrtfnorm,task)
sample = zeros(1,nSample);
i = 1;
if task == 1 % nincs súly, f = lambda*exp(-lambda*x)/fnorm
    while i<= nSample
        preSample = -log(1 - fnorm*rand)/lambda;
        if preSample <= a
            sample(i) = preSample;
            i = i + 1;
        end
    end
elseif task == 2 % nincs súly, f = x/fnorm
    while i<= nSample
        preSample = sqrt(2*rand*fnorm);
        if preSample <= a
            sample(i) = preSample;
            i = i + 1;
        end
    end
elseif task == 3 % p = 1/a, f = lambda*exp(-lambda*x)
    while i<= nSample
        preSample = a*rand;
        if preSample <= a
            sample(i) = preSample;
            i = i + 1;
        end
    end
elseif task == 4 % p = f*g/integral, f = lambda*exp(-lambda*x)
    while i<= nSample
        preSample1 = -log(1 - sqrt(lambda*fnorm)*rand*sqrtfnorm)/lambda;
        preSample2 = -log(1 - sqrt(lambda*fnorm)*rand*sqrtfnorm)/lambda;
        preSample = preSample1 + preSample2;
        if preSample <= a
            sample(i) = preSample;
            i = i + 1;
        end
    end
end
end

%normálásos eljárás (0)///////////////////////////////////////////////////////
function [sample] = calculateSampleNorm(lambda,a,nSample,fnorm,sqrtfnorm,task)
if task == 1 % nincs súly, f = lambda*exp(-lambda*x)/fnorm
    for i = 1:nSample
        sample(i) = -log(1 - fnorm*rand)/lambda;
    end
elseif task == 2 % nincs súly, f = x/fnorm
    for i = 1:nSample
        sample(i) = sqrt(rand*2*fnorm);
    end
elseif task == 3 % p = 1/a, f = lambda*exp(-lambda*x)
    for i = 1:nSample
        sample(i) = a*rand*fnorm;
    end
elseif task == 4 %p = f*g/integral, f = lambda*exp(-lambda*x) (rossz)
    for i = 1:nSample
        sample1 = -log(1 - sqrt(lambda*fnorm)*rand*sqrtfnorm)/lambda;
        sample2 = -log(1 - sqrt(lambda*fnorm)*rand*sqrtfnorm)/lambda;
        sample(i) = sample1 + sample2;
       
    end
end
end
%mintavételezés P bõl ha van súly


