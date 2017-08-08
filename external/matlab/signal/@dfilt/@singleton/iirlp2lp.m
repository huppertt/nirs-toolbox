function [Ht, anum, aden] = iirlp2lp(Ho, varargin)
%IIRLP2LP IIR lowpass to lowpass frequency transformation.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[Ht, anum, aden] = iirxform(Ho, @iirlp2lp, varargin{:});

% [EOF]
