function [F,A,W] = super_getNumericSpecs(h,d)
%GETNUMERICSPECS  Get and evaluate design specs from object.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.



F = get(d,'FrequencyVector');
A = get(d,'MagnitudeVector');
W = get(d,'WeightVector');

