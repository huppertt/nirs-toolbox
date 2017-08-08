function [Ht, anum, aden] = iirlp2mb(Hd, varargin)
%IIRLP2MB IIR Lowpass to multiband transformation

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[Ht, anum, aden] = iirxform(Hd, @iirlp2mb, varargin{:});

% [EOF]
