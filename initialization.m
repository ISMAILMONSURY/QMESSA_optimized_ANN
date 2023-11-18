
%% initialization

function [x,fit]=initialization(pop,dim,ub,lb,fobj)

xxx=rand(1);
YY=ones(1,dim);
for i =1:dim
    xxx=mod(xxx+0.4204-(0.0305/pi)*sin(2*pi*xxx),1);
    YY(i)=xxx;   
end

x=ones(pop,dim);
fit=ones(1,pop);
x2=x;
fit2=fit;
byx1=x;
byx2=x;

for i = 1 : pop
    
    x(i,:) = lb + (ub - lb) .* YY;   
    fit(i) = fobj(x(i,:)) ;  
    theta=asin(2*(x(i,:)-lb)./(ub-lb)-1);  
    x2(i,:)=0.5*(1+cos(theta)).*(ub-lb)+lb;
    fit2(i)=fobj(x2(i,:));
    if fit2(i)<fit(i)   
        x(i,:)=x2(i,:);
        fit(i)=fit2(i);
    end
    theta00=rand*(pi/4)+(pi/4);
    H=[cos(theta00) -sin(theta00);sin(theta00) cos(theta00)];
    COS=cos(theta);
    SIN=sin(theta);
    bit=[COS;SIN];
    bitby=H*bit;  
    byx1(i,:)=0.5*(1+bitby(1,:)).*(ub-lb)+lb; 
    byx2(i,:)=0.5*(1+bitby(2,:)).*(ub-lb)+lb; 
    byx1(i,:)=Bounds(byx1(i,:),lb,ub);
    byx2(i,:)=Bounds(byx2(i,:),lb,ub);
    fit21=fobj(byx1(i,:));
    fit22=fobj(byx2(i,:));
    if fit21<fit(i) && fit22>=fit(i)
        x(i,:)=byx1(i,:);
        fit(i)=fit21;
    elseif fit22<fit(i) && fit21>=fit(i)
        x(i,:)=byx2(i,:);
        fit(i)=fit22;
    elseif fit21<fit(i) && fit22<fit(i)
        if fit21<fit22
            x(i,:)=byx1(i,:);
            fit(i)=fit21;
        else
            x(i,:)=byx2(i,:);
            fit(i)=fit22;
        end
    end  
end

function s = Bounds(s, Lb, Ub)
  temp = s;
  I = temp < Lb;     
  temp(I) = Lb(I);   
  J = temp > Ub;       
  temp(J) = Ub(J);  
  s = temp; 
end

end