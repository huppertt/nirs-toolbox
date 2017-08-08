function [p, v] = coefficient_info(this)
%COEFFICIENT_INFO   

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

p = {getString(message('signal:dfilt:info:NumberofMultipliers'))};
v = {sprintf('%d', length(this.AllpassCoefficients))};

% [EOF]

