function [s,str] = createStruct(h)
%CREATESTRUCT Create a structure for each type.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Return a string for the name of the enumerated data type that will be created
str = 'butterFilterTypes';

s(1).construct.specify = 'filtdes.lpbutter';
s(1).construct.minimum = 'filtdes.lpminbutter';
s(1).tag = 'Lowpass';

s(2).construct.specify = 'filtdes.hpbutter';
s(2).construct.minimum = 'filtdes.hpminbutter';
s(2).tag = 'Highpass';

s(3).construct.specify = 'filtdes.bpbutter';
s(3).construct.minimum = 'filtdes.bpminbutter';
s(3).tag = 'Bandpass';

s(4).construct.specify = 'filtdes.bsbutter';
s(4).construct.minimum = 'filtdes.bsminbutter';
s(4).tag = 'Bandstop';
