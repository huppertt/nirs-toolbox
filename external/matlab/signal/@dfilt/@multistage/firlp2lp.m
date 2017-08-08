function Ht = firlp2lp(Hd, varargin)
%FIRLP2LP FIR Lowpass to lowpass transformation

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

Ht = firxform(Hd, @firlp2lp, varargin{:});

% [EOF]
