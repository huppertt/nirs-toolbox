function [F,E,A,W] = getNumericSpecs(h,d)
%GETNUMERICSPECS  Get and evaluate design specs from object.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.


Fpass1 = get(d,'Fpass1');
Fstop1 = get(d,'Fstop1');
Fstop2 = get(d,'Fstop2');
Fpass2 = get(d,'Fpass2');

Wpass1 = get(d,'Wpass1');
Wstop = get(d,'Wstop');
Wpass2 = get(d,'Wpass2');

Fedges = [Fpass1, Fstop1, Fstop2, Fpass2];

% Form vectors of args
F = [0, Fedges, 1];
E = F;
A = [1 1 0 0 1 1];
W = [Wpass1, Wpass1, Wstop, Wstop, Wpass2, Wpass2];
