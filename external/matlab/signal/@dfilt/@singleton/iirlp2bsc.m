function [Ht, anum, aden] = iirlp2bsc(Ho, varargin)
%IIRLP2BSC IIR lowpass to complex bandstop frequency transformation.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[Ht, anum, aden] = iirxform(Ho, @iirlp2bsc, varargin{:});

% [EOF]
