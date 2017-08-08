function [s,str] = createStruct(h)
%CREATESTRUCT Create a structure for each type.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Return a string for the name of the enumerated data type that will be created
str = 'cheby2FilterTypes';

s(1).construct.specify = 'filtdes.lpcheby2';
s(1).construct.minimum = 'filtdes.lpmincheby2';
s(1).tag = 'Lowpass';

s(2).construct.specify = 'filtdes.hpcheby2';
s(2).construct.minimum = 'filtdes.hpmincheby2';
s(2).tag = 'Highpass';

s(3).construct.specify = 'filtdes.bpcheby2';
s(3).construct.minimum = 'filtdes.bpmincheby2';
s(3).tag = 'Bandpass';

s(4).construct.specify = 'filtdes.bscheby2';
s(4).construct.minimum = 'filtdes.bsmincheby2';
s(4).tag = 'Bandstop';