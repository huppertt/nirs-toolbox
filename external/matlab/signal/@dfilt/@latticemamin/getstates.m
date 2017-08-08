function S = getstates(Hm,S)
%GETSTATES Overloaded get for the States property.

% This should be a private method

%   Author: R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.

S = Hm.HiddenStates;
S = S(1:end-2,:);

% Make sure we return a true [].
if isempty(S)
    S = [];
end

% [EOF]
