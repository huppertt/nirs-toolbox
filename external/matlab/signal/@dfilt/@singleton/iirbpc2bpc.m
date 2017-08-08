function [Ht, anum, aden] = iirbpc2bpc(Ho, varargin)
%IIRBPC2BPC IIR complex bandpass to complex bandpass frequency transformation.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[Ht, anum, aden] = iirxform(Ho, @iirbpc2bpc, varargin{:});

% [EOF]
