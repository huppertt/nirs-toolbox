function [F,A,W] = getNumericSpecs(h,d)
%GETNUMERICSPECS  Get and evaluate design specs from object.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Call alternate method
[F,A,W] = super_getNumericSpecs(h,d);
