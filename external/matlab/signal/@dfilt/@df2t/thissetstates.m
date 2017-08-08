function S = thissetstates(Hd,S)
%THISSETSTATES Overloaded set for the States property.

% This should be a private method

%   Author: V. Pellissier
%   Copyright 1988-2003 The MathWorks, Inc.

% No longer need to append zeros to states according to new df2t
% filter implementation

% Check data type, quantize if needed
S = validatestates(Hd.filterquantizer, S);

Hd.HiddenStates = S;
