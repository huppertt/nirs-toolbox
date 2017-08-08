function S = thissetstates(Hd,S)
%THISSETSTATES Overloaded set for the States property.

% This should be a private method

%   Author: R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

if isempty(S),
    S = nullstate1(Hd.filterquantizer);
else
    % Check data type, quantize if needed
    S = validatestates(Hd.filterquantizer, S);
    S = [S;zeros(1,size(S,2))];
end

Hd.HiddenStates = S;
