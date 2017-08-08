function [s,str] = createStruct(h)
%CREATESTRUCT Create a structure for each type.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

% Return a string for the name of the enumerated data type that will be created
str = 'gremezFilterTypes';

% Add the minimum order version for each of the response types

s(1).construct.specify = 'filtdes.gremezlp';
s(1).construct.minimum = 'filtdes.gremezlpmin';
s(1).tag = 'Lowpass';

s(2).construct.specify = 'filtdes.gremezhp';
s(2).construct.minimum = 'filtdes.gremezhpmin';
s(2).tag = 'Highpass';

s(3).construct.specify = 'filtdes.gremezbp';
s(3).construct.minimum = 'filtdes.gremezbpmin';
s(3).tag = 'Bandpass';

s(4).construct.specify = 'filtdes.gremezbs';
s(4).construct.minimum = 'filtdes.gremezbsmin';
s(4).tag = 'Bandstop';

s(5).construct.specify = 'filtdes.gremezarbmag';
s(5).construct.minimum = 'filtdes.gremezarbmagmin';
s(5).tag = 'Arbitrary magnitude';

s(6).construct.specify = 'filtdes.gremezdiff';
s(6).construct.minimum = '';
s(6).tag = 'Differentiator';

s(7).construct.specify = 'filtdes.gremezhilb';
s(7).construct.minimum = '';
s(7).tag = 'Hilbert transformer';

s(8).construct.specify = 'filtdes.gremezmultiband';
s(8).construct.minimum = 'filtdes.gremezmultibandmin';
s(8).tag = 'Multiband';

% [EOF]
