function e = classiferror(C,Sfit,W,~)

%   Copyright 2010-2013 The MathWorks, Inc.

% If one class, return zero error
if size(C,2)==1
    e = 0;
    return;
end

% Classification error is the fraction of incorrect predictions
notNaN = ~all(isnan(Sfit),2);
[~,y] = max(C(notNaN,:),[],2);
[~,yfit] = max(Sfit(notNaN,:),[],2);
W = W(notNaN,:);
e = sum((y~=yfit).*W) / sum(W);
end
