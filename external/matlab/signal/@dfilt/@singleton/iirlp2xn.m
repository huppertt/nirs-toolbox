function [Ht, anum, aden] = iirlp2xn(Ho, varargin)
%IIRLP2XN IIR lowpass to N-point frequency transformation.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[Ht, anum, aden] = iirxform(Ho, @iirlp2xn, varargin{:});

% [EOF]
