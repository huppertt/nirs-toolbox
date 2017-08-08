function [Ht, anum, aden] = iirlp2hp(Ho, varargin)
%IIRLP2HP IIR lowpass to highpass frequency transformation.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[Ht, anum, aden] = iirxform(Ho, @iirlp2hp, varargin{:});

% [EOF]
