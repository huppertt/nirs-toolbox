function [s,str] = createStruct(h)
%CREATESTRUCT Create a structure for each type.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

% Return a string for the name of the enumerated data type that will be created
str = 'fircbandFilterTypes';

s(1).construct.specify = 'filtdes.fircbandlp';
s(1).construct.c1 = 'filtdes.fircbandlppass';
s(1).construct.c2 = 'filtdes.fircbandlpstop';
s(1).tag = 'Lowpass';

s(2).construct.specify = 'filtdes.fircbandhp';
s(2).construct.c1 = 'filtdes.fircbandhpstop';
s(2).construct.c2 = 'filtdes.fircbandhppass';
s(2).tag = 'Highpass';

s(3).construct.specify = 'filtdes.fircbandbp';
s(3).construct.c1  = 'filtdes.fircbandbpstop1';
s(3).construct.c2  = 'filtdes.fircbandbppass';
s(3).construct.c3  = 'filtdes.fircbandbpstop2';
s(3).construct.c12 = 'filtdes.fircbandbpstoppass';
s(3).construct.c13 = 'filtdes.fircbandbpstop';
s(3).construct.c23 = 'filtdes.fircbandbppassstop';
s(3).tag = 'Bandpass';

s(4).construct.specify = 'filtdes.fircbandbs';
s(4).construct.c1  = 'filtdes.fircbandbspass1';
s(4).construct.c2  = 'filtdes.fircbandbsstop';
s(4).construct.c3  = 'filtdes.fircbandbspass2';
s(4).construct.c12 = 'filtdes.fircbandbspassstop';
s(4).construct.c13 = 'filtdes.fircbandbspass';
s(4).construct.c23 = 'filtdes.fircbandbsstoppass';
s(4).tag = 'Bandstop';

s(5).construct.specify = 'filtdes.fircbandarbmag';
s(5).tag = 'Arbitrary Magnitude';

s(6).construct.specify = 'filtdes.fircbanddiff';
s(6).tag = 'Differentiator';

s(7).construct.specify = 'filtdes.fircbandmultiband';
s(7).tag = 'Multiband';

% [EOF]
