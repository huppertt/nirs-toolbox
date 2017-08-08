function [Ht,anum,aden] = iirlp2bsc(Hd, varargin)
%IIRLP2BSC IIR Lowpass to complex bandstop transformation

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.

[Ht,anum,aden] = ciirxform(Hd, @zpklp2bsc, varargin{:});

