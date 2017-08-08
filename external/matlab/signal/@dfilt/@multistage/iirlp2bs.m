function [Ht, anum, aden] = iirlp2bs(Hd, varargin)
%IIRLP2BS IIR Lowpass to bandstop transformation

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[Ht, anum, aden] = iirxform(Hd, @iirlp2bs, varargin{:});

% [EOF]
