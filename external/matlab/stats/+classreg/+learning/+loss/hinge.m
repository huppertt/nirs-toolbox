function loss = hinge(C,Sfit,W,cost)

%   Copyright 2013-2013 The MathWorks, Inc.

[~,y] = max(C,[],2);
[N,K] = size(C);
Sfit = Sfit(sub2ind([N K],(1:N)',y));
notNaN = ~isnan(Sfit);

hloss = 1 - Sfit(notNaN);
hloss = max(0,hloss);

loss = sum( W(notNaN).*hloss ) / sum(W(notNaN));
end
