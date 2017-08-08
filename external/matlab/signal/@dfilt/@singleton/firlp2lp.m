function Ht = firlp2lp(Ho)
%FIRLP2LP FIR Lowpass to lowpass frequency transformation

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

if ~isfir(Ho) | ~islinphase(Ho) | firtype(Ho) ~= 1,
    error(message('signal:dfilt:singleton:firlp2lp:DFILTErr'));
end

% Special case for scalars
if isscalar(Ho),
    % Ensure that the gain is preserved 
    Ht = copy(Ho);
    return;
end
    
Ht = firxform(Ho, @firlp2lp);

% [EOF]
