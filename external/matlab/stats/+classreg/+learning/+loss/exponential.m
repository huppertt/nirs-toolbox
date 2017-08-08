function loss = exponential(C,Sfit,W,cost)

%   Copyright 2010 The MathWorks, Inc.



[~,y] = max(C,[],2);
[N,K] = size(C);
Sfit = Sfit(sub2ind([N K],(1:N)',y));
notNaN = ~isnan(Sfit);
loss = sum( W(notNaN).*exp(-Sfit(notNaN)) ) / sum(W(notNaN));
end
