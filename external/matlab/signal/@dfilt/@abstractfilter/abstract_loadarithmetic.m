function abstract_loadarithmetic(this, s)
%ABSTRACT_LOADARITHMETIC   Load the arithmetic information.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

% If s is not a structure we need to copy the privfq.
if ~isstruct(s)
    
    for indx = 1:length(s.privfq)
        privfq(indx) = copy(s.privfq(indx));
    end
    arith = s.privArithmetic;
else
    
    if isfield(s, 'privArithmetic')
        arith = s.privArithmetic;
    else
        arith = s.Arithmetic;
    end
    privfq = s.privfq;
end

set(this, 'privfq', privfq, ...
    'privArithmetic', arith);

% [EOF]
