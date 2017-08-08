function [p, v] = coefficient_info(this)
%COEFFICIENT_INFO   

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

p = {getString(message('signal:dfilt:info:NumberofSections'))};
v = {sprintf('%d', nsections(this))};

% [EOF]
