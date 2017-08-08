function [s,str] = createStruct(h)
%CREATESTRUCT Create a structure for each type.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Return a string for the name of the enumerated data type that will be created
str = 'cheby1FilterTypes';

s(1).construct.specify = 'filtdes.lpcheby1';
s(1).construct.minimum = 'filtdes.lpmincheby1';
s(1).tag = 'Lowpass';

s(2).construct.specify = 'filtdes.hpcheby1';
s(2).construct.minimum = 'filtdes.hpmincheby1';
s(2).tag = 'Highpass';

s(3).construct.specify = 'filtdes.bpcheby1';
s(3).construct.minimum = 'filtdes.bpmincheby1';
s(3).tag = 'Bandpass';

s(4).construct.specify = 'filtdes.bscheby1';
s(4).construct.minimum = 'filtdes.bsmincheby1';
s(4).tag = 'Bandstop';

