function [Ht,anum,aden] = iirlp2bpc(Hd, varargin)
%IIRLP2BPC IIR Lowpass to complex bandpass transformation

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.

[Ht,anum,aden] = ciirxform(Hd, @zpklp2bpc, varargin{:});
