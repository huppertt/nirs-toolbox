function [yfilt,f] = innovations(Y,Pmax,verbose,maxBICSearch)
% This removes auto-correlation and returns the innvations model;

[m, n] = size(Y);
Pmax = min(m,round(Pmax));

if(nargin<3)
    verbose=false;
end
if nargin < 4 || isempty(maxBICSearch), maxBICSearch = 0; end

% Cap effective Pmax for BIC search
if maxBICSearch > 0
    effectivePmax = min(Pmax, maxBICSearch);
else
    effectivePmax = Pmax;
end

% model selection
yfilt = nan(size(Y));
f = cell(1,n);

if verbose
    h=waitbar(0,'Computing AR-model');
    for i = 1:n
        h=waitbar(i/n,h);
        y1 = mean(Y(:,i));
        y = bsxfun(@minus,Y(:,i),y1);
        a = nirs.math.ar_fit(y, effectivePmax);
        f{i}=[1; -a(2:end)];
        yfilt(:,i) = filter(f{i}, 1, y);
    end
    try; close(h); end;
else
    parfor i = 1:n
        y1 = mean(Y(:,i));
        y = bsxfun(@minus,Y(:,i),y1);
        a = nirs.math.ar_fit(y, effectivePmax);
        f{i}=[1; -a(2:end)];
        yfilt(:,i) = filter(f{i}, 1, y);
    end
end

return