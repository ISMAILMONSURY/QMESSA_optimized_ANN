% QMESSA
function [bestX,fMin] = QMESSA(dim,fobj)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializing the parameters of improved sparrow search algorithm %

pop = 50;    % searchagents_no
M   =300;    % maxiteration
c   = -2;     %lb (lower bound)
d   = 2;      %ub (upper bound)
lb  = -2;
ub  = 2;


ST=0.8;
P_percent = 0.2;
S_percent = 0.2;
PD = round(pop*P_percent);
SD = round(pop*S_percent);
lb = c.*ones(1,dim);
ub = d.*ones(1,dim);
m=M/5;
a=exp(log(m)/(m));

% initialization
[x,fit]=initialization(pop,dim,ub,lb,fobj);

pFit = fit;
pX = x;       % The individual's best position corresponding to the pFit
[fMin, bestI] = min(fit);
bestX = x(bestI,:);

% Start updating the solutions.
for t = 1 : M
    
    [~, sortIndex] = sort(pFit);
    [fMin, bestI] = min(pFit);
    bestX = pX(bestI,:);
    bestX2=pX(sortIndex(2),:);
    [fMax,worseI]=max(pFit);
    worseX= pX(worseI,:);
    r2=rand(1);
    
    RR=1-t/M;
    
    Q1=trnd(a^t/10,[1,PD]);
    for i = 1 : PD
        if(r2 < ST)
            alpha=rand(1);
            miu=rand*1+1.5;tao=rand*8+4;
            x(sortIndex(i),:) = pX(sortIndex(i),:)*(    miu/(exp((i^tao)/alpha/M))*(2*rand-1)*rand^0.5+(rand-0.5)^3  );
            x(sortIndex(i),:) = Bound(x(sortIndex(i),:),lb,ub,bestX,bestX2);
            fit(sortIndex(i)) = fobj(x(sortIndex(i),:));
        else
            
            L=ones(1,dim);
            x(sortIndex(i),:) = pX(sortIndex(i),:)+Q1(i)*L;
            x(sortIndex(i),:) = Bound(x(sortIndex(i),:),lb,ub,bestX,bestX2);
            fit(sortIndex(i)) = fobj(x(sortIndex(i),:));
        end
    end
    
    [~,bestII] = min(fit);
    bestP = x(bestII,:);
    
    
    Q1=trnd(a^t/10,[1,pop-PD]);
    for i = (PD + 1) : pop
        A=floor(rand(1,dim)*2)*2-1;
        if(i>(pop/2))
            x(sortIndex(i),:)=Q1(i-PD)*exp((worseX-pX(sortIndex(i),:))/(i)^2);
        else
            L=ones(1,dim);
            x(sortIndex(i),:)=bestP+(abs((pX(sortIndex(i),:)-bestP)))*(A'*(A*A')^(-1))*L;
        end
        x(sortIndex(i),:) = Bound(x(sortIndex(i),:),lb,ub,bestX,bestX2);
        fit(sortIndex(i)) = fobj(x(sortIndex(i),:));
    end
    
    c=randperm(numel(sortIndex));
    b=sortIndex(c(1:SD));
    Q1=trnd(a^t/10,[length(b)+1,dim]);
    for j = 1:length(b)
        if(pFit(sortIndex(b(j))) == (fMin))
            K=(2*rand(1)-1);
            epsilon=1e-100;
            x(sortIndex(b(j)),:) =pX(sortIndex(b(j)),:)+K*(abs(pX(sortIndex(b(j)),:)-worseX))/ (pFit(sortIndex(b(j)))-fMax+epsilon);
        else
            x(sortIndex(b(j)),:)=bestX+Q1(j,:).*(abs((pX(sortIndex(b(j)),:) - bestX)));
        end
        x(sortIndex(b(j)),:) = Bound(x(sortIndex(b(j)),:),lb,ub,bestX,bestX2);
        fit(sortIndex(b(j))) = fobj(x(sortIndex(b(j)),:));
        
    end
    
    
    for i = 1 : pop
        if (fit(i) < pFit(i))
            pFit(i) = fit(i);
            pX(i,:) = x(i,:);
        end
        
        if(pFit(i) < fMin)
            fMin= pFit(i);
            bestX = pX(i,:);
        end
    end
    
    [~, sortIndex] = sort(pFit);
    Q1=trnd(a^t/10,[1,pop-round(0.8*pop)+1]);
    for ii = round(0.8*pop) : pop
        pX(sortIndex(ii),:)=(ub-lb)*RR*Q1(ii-(round(0.8*pop-1)))+bestX;
        pX(sortIndex(ii),:) = Bound(pX(sortIndex(ii),:),lb,ub,bestX,bestX2);
    end

    Convergence_curve(t)=fMin;
    
end

end

