function [yfilt,f] = innovations(Y,Pmax,verbose)
% This removes auto-correlation and returns the innvations model;

[m, n] = size(Y);
Pmax = min(m,Pmax);
X = []; lst = [];

if(nargin<3)
    verbose=false;
end

if(verbose)
    h=waitbar(0,'Computing AR-model');
end


% model selection
yfilt = nan(size(Y));
f = cell(1,size(Y,2));
for i = 1:n
    % fit AR model
    
    
    if(verbose)
        h=waitbar(i/n,h);
    end
    
    y1 = mean(Y(:,i));
    y = bsxfun(@minus,Y(:,i),y1);
    a = nirs.math.ar_fit(y, Pmax);
    f{i}=[1; -a(2:end)];
    
    yfilt(:,i) = filter(f{i}, 1, y);
    
end


if(verbose)
    try; close(h); end;
end

return