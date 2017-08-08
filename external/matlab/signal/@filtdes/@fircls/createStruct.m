function [s,str] = createStruct(h)
%CREATESTRUCT Create a structure for each type.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

% Return a string for the name of the enumerated data type that will be created
str = 'firclsFilterTypes';

s(1).construct = 'filtdes.firclslp';
s(1).tag = 'Lowpass';

s(2).construct = 'filtdes.firclshp';
s(2).tag = 'Highpass';

s(3).construct = 'filtdes.firclsbp';
s(3).tag = 'Bandpass';

s(4).construct = 'filtdes.firclsbs';
s(4).tag = 'Bandstop';

s(5).construct = 'filtdes.firclsmultiband';
s(5).tag       = 'Multiband';

% [EOF]
