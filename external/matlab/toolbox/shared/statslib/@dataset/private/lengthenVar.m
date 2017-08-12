function b = lengthenVar(a,n)
% Lengthen an existing variable out to n rows.

%   Copyright 2006-2013 The MathWorks, Inc.

m = size(a,1);
b = a;
if isnumeric(a)
    b(m+1:n,:) = 0;
elseif islogical(a)
    b(m+1:n,:) = false;
elseif isa(a,'categorical')
    b(m+1:n,:) = categorical.undefLabel;
elseif iscell(a)
    b(m+1:n,:) = {[]};
else % including struct and objects
    if ismatrix(a)
        b(n+1,:) = b(1,:); b = b(1:n,:); % without using reshape, may not be one
    else
        sizeOut = size(a); sizeOut(1) = n;
        b(n+1,:) = b(1,:); b = reshape(b(1:n,:),sizeOut);
    end
end
