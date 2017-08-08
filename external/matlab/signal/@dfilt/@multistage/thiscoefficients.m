function c = thiscoefficients(Hd)
%THISCOEFFICIENTS Filter coefficients.
%   C = THISCOEFFICIENTS(Hd) returns a cell array of coefficients of
%   discrete-time filter Hd.
%
%   See also DFILT.   

%   Author: Thomas A. Bryan
%   Copyright 1988-2004 The MathWorks, Inc.

c = {};
for k=1:length(Hd.Stage)
    c = {c{:}, coefficients(Hd.Stage(k))};
end
