function [F,E,A,W] = getNumericSpecs(h,d)
%GETNUMERICSPECS  Get and evaluate design specs from object.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

[F,A,W] = super_getNumericSpecs(h,d);

E = get(d,'FrequencyEdges');
