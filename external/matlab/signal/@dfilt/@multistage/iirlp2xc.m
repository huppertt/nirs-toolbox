function [Ht, anum, aden] = iirlp2xc(Hd, varargin)
%IIRLP2XC IIR Lowpass to complex N-Point transformation

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[Ht, anum, aden] = iirxform(Hd, @iirlp2xc, varargin{:});

% [EOF]
