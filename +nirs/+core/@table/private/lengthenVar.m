function b = lengthenVar(a,n)
% LENGTHENVAR Lengthen an existing variable out to n rows.

%   Copyright 2013-2014 The MathWorks, Inc.

m = size(a,1);
b = a;
if isnumeric(a)
    % let a numeric subclass pad with its choice, e.g. NaN
    b(n+1,:) = 0; b(n+1,:) = [];
elseif islogical(a)
    b(n,:) = false;
elseif isa(a,'categorical')
    b(n,:) = categorical.undefLabel;
elseif iscell(a)
    b(m+1:n,:) = {[]};
elseif istable(a)
    % can't use table subscripting directly
    b = subsasgnParens(b,{n+1 ':'},subsrefParens(b,{1 ':'}));
    b = subsrefParens(b,{1:n ':'});
else % including struct and objects
    if ismatrix(a)
        b(n+1,:) = b(1,:); b = b(1:n,:); % without using reshape, may not be one
    else
        sizeOut = size(a); sizeOut(1) = n;
        b(n+1,:) = b(1,:); b = reshape(b(1:n,:),sizeOut);
    end
end
