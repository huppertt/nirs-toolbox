function [s,str] = createStruct(h)
%CREATESTRUCT Create a structure for each type.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Return a string for the name of the enumerated data type that will be created
str = 'ellipFilterTypes';

s(1).construct.specify = 'filtdes.lpellip';
s(1).construct.minimum = 'filtdes.lpminellip';
s(1).tag = 'Lowpass';

s(2).construct.specify = 'filtdes.hpellip';
s(2).construct.minimum = 'filtdes.hpminellip';
s(2).tag = 'Highpass';

s(3).construct.specify = 'filtdes.bpellip';
s(3).construct.minimum = 'filtdes.bpminellip';
s(3).tag = 'Bandpass';

s(4).construct.specify = 'filtdes.bsellip';
s(4).construct.minimum = 'filtdes.bsminellip';
s(4).tag = 'Bandstop';
