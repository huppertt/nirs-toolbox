function [s,str] = createStruct(h)
%CREATESTRUCT Create a structure for each type.

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

% Return a string for the name of the enumerated data type that will be created
str = 'cremezFilterTypes';

s(1).construct = 'filtdes.cremezlp';
s(1).tag       = 'Lowpass';

s(2).construct = 'filtdes.cremezhp';
s(2).tag       = 'Highpass';

s(3).construct = 'filtdes.cremezbp';
s(3).tag       = 'Bandpass';

s(4).construct = 'filtdes.cremezbs';
s(4).tag       = 'Bandstop';

s(5).construct = 'filtdes.cremezmultiband';
s(5).tag       = 'Multiband';

s(6).construct = 'filtdes.cremezhilb';
s(6).tag       = 'Hilbert Transformer';

s(7).construct = 'filtdes.cremezdiff';
s(7).tag       = 'Differentiator';
 
s(8).construct = 'filtdes.cremezlpinvsinc';
s(8).tag       = 'Inverse Sinc Lowpass';

% [EOF]
