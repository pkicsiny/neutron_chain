function [x, f] = calculatePDF2(sample, nPoints, normalize)

%intervallumok
s=(max(sample)-min(sample))/nPoints;
for i=1:nPoints
    x(i) = min(sample)+(2*i-1)/2*s; 
end
%minta elemeinek berakása az intervallumokba
f(1)=0;
for j=1:length(sample)
    if sample(j)>=min(sample) && sample(j)<=min(sample)+s;
        f(1)=f(1)+1;
    end
end

for i=2:nPoints
    f(i)=0;
    for j=1:length(sample)
        if sample(j)>x(i)-s/2 && sample(j)<=x(i)+s/2
            f(i)=f(i)+1;
        end
    end
    
end
if normalize==1
   
        f=f/(sum(f)*s);
end
end

