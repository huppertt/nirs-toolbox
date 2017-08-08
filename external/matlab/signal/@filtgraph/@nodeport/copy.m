function NP = copy(np)
%COPY Copy method for Nodeport
%   Copy method to force a deep copy of a nodeport vector.

%   Author(s): S Dhoorjaty
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(1,1,nargin,'struct'));
for k = 1:size(np,1)
    for l = 1:size(np,2)
        NP(k,l) = feval(str2func(class(np)));
        NP(k,l).node = np(k,l).node;
        NP(k,l).port = np(k,l).port;
    end
end

