function [s,str] = createStruct(h)
%CREATESTRUCT Create a structure for each type.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Return a string for the name of the enumerated data type that will be created
str = 'maxflatFilterTypes';

s(1).construct = 'filtdes.lpiirmaxflat';
s(1).tag = 'Lowpass';


