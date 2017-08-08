function [s,str] = createStruct(h)
%CREATESTRUCT Create a structure for each type.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Return a string for the name of the enumerated data type that will be created
str = 'lpnormFilterTypes';

s(1).construct = 'filtdes.lpnormlppassstop';
s(1).tag = 'Lowpass';

s(2).construct = 'filtdes.lpnormhppassstop';
s(2).tag = 'Highpass';

s(3).construct = 'filtdes.lpnormbppassstop';
s(3).tag = 'Bandpass';

s(4).construct = 'filtdes.lpnormbspassstop';
s(4).tag = 'Bandstop';

s(5).construct = 'filtdes.freqedgemagweight';
s(5).tag = 'Arbitrary magnitude';