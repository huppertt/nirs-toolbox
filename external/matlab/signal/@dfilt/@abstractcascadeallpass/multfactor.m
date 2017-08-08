function [f offset] = multfactor(this)
%MULTFACTOR   

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.


f = ones(1,length(this.refallpasscoeffs));
offset = zeros(size(f));

% [EOF]
