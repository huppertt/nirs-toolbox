function [yfilt,f] = innovations(Y,Pmax)
% This removes auto-correlation and returns the innvations model;

[m, n] = size(Y);

X = []; lst = [];

% model selection
for i = 1:n
    % fit AR model
    y1 = mean(Y(:,1));
    y = bsxfun(@minus,Y(:,i),y1);
    a = nirs.math.ar_fit(y, Pmax);
    f{i}=[1; -a(2:end)];
    
    yfilt(:,i) = filter(f{i}, 1, y);
    
end

return