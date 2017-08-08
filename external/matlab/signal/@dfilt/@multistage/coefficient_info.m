function [p,v] = coefficient_info(this)
%COEFFICIENT_INFO   

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

p = {getString(message('signal:dfilt:info:NumberofStages'))};
v = {num2str(nstages(this))};

% [EOF]
