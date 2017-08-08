function [Ht,anum,aden] = iirlp2xc(Hd, varargin)
%IIRLP2XC IIR Lowpass to complex N-Point transformation

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.

[Ht,anum,aden] = ciirxform(Hd, @zpklp2xc, varargin{:});

