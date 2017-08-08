function [yfilt, ARIMA_mdl] = innovations_ARIMA(Y,Pmax,nMA,nI,verbose)
% This removes auto-correlation and returns the innvations model;

[m, n] = size(Y);

X = []; lst = [];


if(nargin<4)
    nI=0;
end

if(nargin<5)
    verbose=false;
end

if(verbose)
    h=waitbar(0,'Computing AR-model');
end


% model selection
yfilt = Y;
ARIMA_mdl = cell(1,size(Y,2));
for i = 1:n
    % fit AR model
    
    if(verbose)
        h=waitbar(i/n,h);
    end
    
    y1 = mean(Y(:,i));
    y = bsxfun(@minus,Y(:,i),y1);
    
    mdl=arima(Pmax,nI,nMA);
    mdl=estimate(mdl,y,'Display','off');
    yfilt(:,i)=mdl.infer(y);
    ARIMA_mdl{i}=mdl;
end


if(verbose)
    try; close(h); end;
end

return