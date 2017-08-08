function validatercf(this, ifactor, dfactor)
%VALIDATERCF   Validate the rcf

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

if nargin < 3
    dfactor = get(this, 'DecimationFactor');
    if nargin < 2
        ifactor = get(this, 'InterpolationFactor');
    end
end

g = gcd(ifactor,dfactor);

if g > 1
    error(message('signal:fdesign:rsrc:validatercf:filterdesign:srconvert:primefactors', 'L', 'M', 'InterpolationFactor', ifactor/g, 'DecimationFactor', dfactor/g)); 
end


% [EOF]
