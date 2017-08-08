function flag = mustNormalize(h)
%MUSTNORMALIZE Method to determine if frequencies must be normalized.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

flag = 0;

freqUnitsOpts = set(h,'freqUnits');
if ~strcmpi(freqUnitsOpts{1},get(h,'freqUnits')), % Not normalized
    flag = 1;
end

