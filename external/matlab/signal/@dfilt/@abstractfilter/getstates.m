function S = getstates(Hm,S)
%GETSTATES Overloaded get for the States property.

% This should be a private method

%   Author: V. Pellissier
%   Copyright 1988-2002 The MathWorks, Inc.

S = Hm.HiddenStates;
