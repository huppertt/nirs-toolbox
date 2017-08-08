function C = classCount(expectedY,observedY)
%CLASSCOUNT Match expected and observed class labels.
%   C=CLASSCOUNT(EXPY,OBSY) returns an N-by-K logical matrix for N
%   observations in vector of observed labels OBSY and K values in vector
%   of expected labels EXPY. The values in EXPY must be unique. Element
%   C(I,J) is true if OBSY(I) equals EXPY(J) and false otherwise.
%
%   If OBSY has elements not found in EXPY, CLASSCOUNT throws an error.
%   Note that this will occur if either argument contains a NaN.
    
%   Copyright 2010-2014 The MathWorks, Inc.

K = length(expectedY);
N = length(observedY);
C = false(N,K);

[tf,grp] = ismember(observedY,expectedY);
if ~all(tf)
    idx = find(~tf,1,'first');
    if     isa(observedY,'classreg.learning.internal.ClassLabel') ...
            || isa(observedY,'categorical') || iscellstr(observedY)
        str = char(observedY(idx));
    else
        str = num2str(observedY(idx,:));
    end
    if isa(observedY,'classreg.learning.internal.ClassLabel')
        cls = class(labels(observedY));
    else
        cls = class(observedY);
    end
    error(message('stats:classreg:learning:internal:classCount:UnknownClass', str, cls));
end

C(sub2ind([N K],(1:N)',grp)) = true;
end
