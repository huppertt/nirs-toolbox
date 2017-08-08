function designobj = getdesignobj(this, str)
%GETDESIGNOBJ   Get the designobj.

%   Author(s): P. Costa
%   Copyright 2005 The MathWorks, Inc.



%#function fdfmethod.elliphilbertmin
%#function fdfmethod.eqriphilbmin
%#function fdfmethod.iirlinphasehilbertmin
designobj.ellip      = 'fdfmethod.elliphilbertmin';
designobj.equiripple = 'fdfmethod.eqriphilbmin';
designobj.iirlinphase = 'fdfmethod.iirlinphasehilbertmin';


if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
