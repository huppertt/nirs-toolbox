function [Ht, anum, aden] = iirbpc2bpc(Hd, varargin)
%IIRBPC2BPC IIR complex bandpass to complex bandpass transformation

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[Ht, anum, aden] = iirxform(Hd, @iirbpc2bpc, varargin{:});

% [EOF]
