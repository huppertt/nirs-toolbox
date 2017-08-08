function [s,str] = createStruct(h)
%CREATESTRUCT Create a structure for each type.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Return a string for the name of the enumerated data type that will be created
str = 'firlsFilterTypes';

s(1).construct = 'filtdes.firlslppassstop';
s(1).tag = 'Lowpass';

s(2).construct = 'filtdes.firlshppassstop';
s(2).tag = 'Highpass';

s(3).construct = 'filtdes.firlsbppassstop';
s(3).tag = 'Bandpass';

s(4).construct = 'filtdes.firlsbspassstop';
s(4).tag = 'Bandstop';

s(5).construct = 'filtdes.firlsarbmag';
s(5).tag = 'Arbitrary magnitude';

s(6).construct = 'filtdes.firlsdiff';
s(6).tag = 'Differentiator';

s(7).construct = 'filtdes.firlshilb';
s(7).tag = 'Hilbert transformer';

s(8).construct = 'filtdes.firlsmultiband';
s(8).tag = 'Multiband';
