function str = set_arith(h,str)
%SET_ARITH   SetFunction for the Arithmetic property.

%   Author(s): R. Losada
%   Copyright 1988-2005 The MathWorks, Inc.

if strncmpi(str, 'fixed', length(str)) && ~isfixptinstalled
    error(message('signal:dfilt:abstractfilter:set_arith:invalidArithmetic'));
end

h.privArithmetic   = str;
h.privMeasurements = [];

% [EOF]
