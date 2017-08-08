function S = thissetstates(Hd,S)
%THISSETSTATES Overloaded set for the States property.

% This should be a private method

%   Author: V. Pellissier
%   Copyright 1988-2004 The MathWorks, Inc.

S = validatestatesobj(Hd.filterquantizer, S);
Hd.HiddenStates = S;
