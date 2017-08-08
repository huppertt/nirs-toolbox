function [Ht, anum, aden] = iirlp2bpc(Ho, varargin)
%IIRLP2BPC IIR lowpass to complex bandpass frequency transformation.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[Ht, anum, aden] = iirxform(Ho, @iirlp2bpc, varargin{:});

% [EOF]
