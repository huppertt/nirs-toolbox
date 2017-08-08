function Ht = firlp2hp(Ho, varargin)
%FIRLP2HP FIR Lowpass to highpass frequency transformation

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

if ~isfir(Ho)
    error(message('signal:dfilt:singleton:firlp2hp:MustBeFIR'));
end

if nargin > 1 && (~islinphase(Ho) || firtype(Ho) ~= 1) && strcmpi(varargin{1},'wide'),
    error(message('signal:dfilt:singleton:firlp2hp:MustBeLinearPhase'));
end

Ht = firxform(Ho, @firlp2hp, varargin{:});

% [EOF]
