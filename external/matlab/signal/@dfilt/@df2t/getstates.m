function S = getstates(Hd,S)
%GETSTATES Overloaded get for the States property.

% This should be a private method

%   Author: V. Pellissier
%   Copyright 1988-2004 The MathWorks, Inc.

ncoeffs = Hd.ncoeffs;

% If the ncoeffs property is [], then we are in an init condition.
if isempty(ncoeffs)
    return;
end

nb = ncoeffs(1);
na = ncoeffs(2);

S = Hd.HiddenStates;

% Check for the scalar case, i.e., NB & NA = 1. This avoids the states
% returning an empty, it returns a 0x1 double matrix (similar to R13)
if (nb == 1) & (na == 1),
    S = [feval(class(S),zeros(0,1))];
end
