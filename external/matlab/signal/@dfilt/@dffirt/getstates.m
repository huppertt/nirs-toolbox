function S = getstates(Hm,S)
%GETSTATES Overloaded get for the States property.

% This should be a private method

%   Author: R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

S = Hm.HiddenStates;
