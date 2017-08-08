function Ht = firlp2hp(Hd, varargin)
%FIRLP2HP FIR Lowpass to highpass transformation

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

Ht = firxform(Hd, @firlp2hp, varargin{:});

% [EOF]
