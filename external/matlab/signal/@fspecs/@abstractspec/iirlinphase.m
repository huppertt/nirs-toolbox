function  Hd = iirlinphase(this,varargin)
%IIRLINPHASE   IIR quasi linear phase digital filter design.
%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

Hd = design(this, 'iirlinphase', varargin{:});

% [EOF]
