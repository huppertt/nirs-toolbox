function [Ht, anum, aden] = iirlp2xn(Hd, varargin)
%IIRLP2XN IIR Lowpass to N-Point transformation

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[Ht, anum, aden] = iirxform(Hd, @iirlp2xn, varargin{:});

% [EOF]
