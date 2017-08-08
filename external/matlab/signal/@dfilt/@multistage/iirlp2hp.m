function [Ht, anum, aden] = iirlp2hp(Hd, varargin)
%IIRLP2HP IIR Lowpass to highpass transformation

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[Ht, anum, aden] = iirxform(Hd, @iirlp2hp, varargin{:});

% [EOF]
