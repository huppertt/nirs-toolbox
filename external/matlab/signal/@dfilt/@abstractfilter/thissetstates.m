function S = thissetstates(Hm,S)
%THISSETSTATES Overloaded set for the States property.

% This should be a private method

%   Author: V. Pellissier
%   Copyright 1988-2002 The MathWorks, Inc.


Hm.HiddenStates = S;

