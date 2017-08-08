function [Ht, anum, aden] = iirlp2mbc(Hd, varargin)
%IIRLP2MBC IIR Lowpass to complex multiband transformation

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[Ht, anum, aden] = iirxform(Hd, @iirlp2mbc, varargin{:});

% [EOF]
