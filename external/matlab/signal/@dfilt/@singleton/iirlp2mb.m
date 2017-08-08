function [Ht, anum, aden] = iirlp2mb(Ho, varargin)
%IIRLP2MB IIR lowpass to multiband frequency transformation.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[Ht, anum, aden] = iirxform(Ho, @iirlp2mb, varargin{:});

% [EOF]
