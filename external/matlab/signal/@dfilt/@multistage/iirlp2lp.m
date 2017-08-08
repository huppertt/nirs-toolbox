function [Ht,anum,aden] = iirlp2lp(Hd, varargin)
%IIRLP2LP IIR Lowpass to lowpass transformation

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[Ht,anum,aden] = iirxform(Hd, @iirlp2lp, varargin{:});

% [EOF]
