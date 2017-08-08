function [Ht, anum, aden] = iirlp2bp(Ho, varargin)
%IIRLP2BP IIR lowpass to bandpass frequency transformation.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[Ht, anum, aden] = iirxform(Ho, @iirlp2bp, varargin{:});

% [EOF]
