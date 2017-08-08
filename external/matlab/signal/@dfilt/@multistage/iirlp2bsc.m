function [Ht, anum, aden] = iirlp2bsc(Hd, varargin)
%IIRLP2BSC IIR Lowpass to complex bandstop transformation

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[Ht, anum, aden] = iirxform(Hd, @iirlp2bsc, varargin{:});

% [EOF]
