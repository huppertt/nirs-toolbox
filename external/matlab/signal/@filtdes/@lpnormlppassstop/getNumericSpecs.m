function [F,E,A,W] = getNumericSpecs(h,d)
%GETNUMERICSPECS  Get and evaluate design specs from object.

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

args = getarguments(h, d);

F = args{1};
E = F;
A = args{2};
W = [args{3}' args{3}']';
W = W(:)';

% [EOF]
