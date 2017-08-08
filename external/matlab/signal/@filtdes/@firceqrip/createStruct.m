function [s,str] = createStruct(h)
%CREATESTRUCT Create a structure for each type.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Return a string for the name of the enumerated data type that will be created
str = 'firceqripFilterTypes';

s(1).construct = 'filtdes.firceqriplp';
s(1).tag = 'Lowpass';

s(2).construct = 'filtdes.firceqriphp';
s(2).tag = 'Highpass';

s(3).construct = 'filtdes.firceqriplpinvsinc';
s(3).tag = 'Inverse Sinc Lowpass';

s(4).construct = 'filtdes.firceqriphpinvsinc';
s(4).tag = 'Inverse Sinc Highpass';
