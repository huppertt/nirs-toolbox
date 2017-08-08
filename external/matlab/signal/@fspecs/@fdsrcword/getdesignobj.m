function designobj = getdesignobj(this, str)
%GETDESIGNOBJ   Get the design object.

%   Copyright 2007 The MathWorks, Inc.

%#function fdfmethod.lagrangesrc
designobj.lagrange = 'fdfmethod.lagrangesrc';

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
