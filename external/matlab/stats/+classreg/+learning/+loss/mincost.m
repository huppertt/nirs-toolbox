function e = mincost(C,Sfit,W,cost)

%   Copyright 2011-2013 The MathWorks, Inc.


% If no cost is provided, return regular classification error
if isempty(cost)
    K = size(C,2);
    cost = ones(K)-eye(K);
end

% If only zero costs, return 0
if ~any(cost(:))
    e = 0;
    return;
end

% Compute expected cost. Sfit are assumed to be posterior probabilities.
expcost = Sfit*cost;

% Find true and predicted class labels
notNaN = ~all(isnan(expcost),2);
[~,y] = max(C(notNaN,:),[],2);
[~,yfit] = min(expcost(notNaN,:),[],2);
W = W(notNaN,:);

% Compute the incurred cost
e = sum(cost(sub2ind(size(cost),y,yfit)).*W) / sum(W);

end
