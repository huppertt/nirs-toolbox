function [s,str] = createStruct(h)
%CREATESTRUCT Create a structure for each type.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Return a string for the name of the enumerated data type that will be created
str = 'fir1FilterTypes';

s(1).construct.specify = 'filtdes.lpfir1';
s(1).construct.minimum = 'filtdes.lpminfir1';
s(1).tag = 'Lowpass';

s(2).construct.specify = 'filtdes.hpfir1';
s(2).construct.minimum = 'filtdes.hpminfir1';
s(2).tag = 'Highpass';

s(3).construct.specify = 'filtdes.bpfir1';
s(3).construct.minimum = 'filtdes.bpminfir1';
s(3).tag = 'Bandpass';

s(4).construct.specify = 'filtdes.bsfir1';
s(4).construct.minimum = 'filtdes.bsminfir1';
s(4).tag = 'Bandstop';

s(5).construct.specify = 'filtdes.rcosfir1';
s(5).construct.minimum = '';
s(5).tag = 'Raised-cosine';

s(6).construct.specify = 'filtdes.lphalffir1';
s(6).construct.minimum = 'filtdes.lphalfminfir1';
s(6).tag = 'Halfband lowpass';

s(7).construct.specify = 'filtdes.hphalffir1';
s(7).construct.minimum = 'filtdes.hphalfminfir1';
s(7).tag = 'Halfband highpass';

s(8).construct.specify = 'filtdes.nyqfir1';
s(8).construct.minimum = 'filtdes.nyqminfir1';
s(8).tag = 'Nyquist';
