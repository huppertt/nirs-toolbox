function out = ismax(in)

%   Copyright 2010 The MathWorks, Inc.


[N,K] = size(in);
out = zeros(N,K);
[~,colnum] = max(in,[],2);
for k=1:K
    out(colnum==k,k) = 1;
end
end
