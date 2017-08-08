function [Ht, anum, aden] = iirlp2mbc(Ho, varargin)
%IIRLP2MBC IIR lowpass to complex multiband frequency transformation.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[Ht, anum, aden] = iirxform(Ho, @iirlp2mbc, varargin{:});

% [EOF]
