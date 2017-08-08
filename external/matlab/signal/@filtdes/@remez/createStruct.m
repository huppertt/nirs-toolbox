function [s,str] = createStruct(h)
%CREATESTRUCT Create a structure for each type.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Return a string for the name of the enumerated data type that will be created
str = 'remezFilterTypes';

% First element is specify order constructor, second is minimum order
s(1).construct.specify = 'filtdes.remezlppassstop';
s(1).construct.minimum = 'filtdes.remezlpminpassstop';
s(1).tag = 'Lowpass';

s(2).construct.specify = 'filtdes.remezhppassstop';
s(2).construct.minimum = 'filtdes.remezhpminpassstop';
s(2).tag = 'Highpass';

s(3).construct.specify = 'filtdes.remezbppassstop';
s(3).construct.minimum = 'filtdes.remezbpminpassstop';
s(3).tag = 'Bandpass';

s(4).construct.specify = 'filtdes.remezbspassstop';
s(4).construct.minimum = 'filtdes.remezbsminpassstop';
s(4).tag = 'Bandstop';

s(5).construct.specify = 'filtdes.remezarbmag';
s(5).construct.minimum = '';
s(5).tag = 'Arbitrary magnitude';

s(6).construct.specify = 'filtdes.remezdiff';
s(6).construct.minimum = '';
s(6).tag = 'Differentiator';

s(7).construct.specify = 'filtdes.remezhilb';
s(7).construct.minimum = '';
s(7).tag = 'Hilbert transformer';

s(8).construct.specify = 'filtdes.remezmultiband';
s(8).construct.minimum = '';
s(8).tag = 'Multiband';

s(9).construct.specify = 'filtdes.remezlphalf';
s(9).construct.minimum = 'filtdes.remezlphalfmin';
s(9).tag = 'Halfband lowpass';

s(10).construct.specify = 'filtdes.remezhphalf';
s(10).construct.minimum = 'filtdes.remezhphalfmin';
s(10).tag = 'Halfband highpass';

s(11).construct.specify = 'filtdes.remeznyq';
s(11).construct.minimum = '';
s(11).tag = 'Nyquist';

