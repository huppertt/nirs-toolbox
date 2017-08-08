function designobj = getdesignobj(this, str)
%GETDESIGNOBJ   Get the design object.

%   Copyright 2005-2011 The MathWorks, Inc.

%#function fdfmethod.lagrange
designobj.lagrange = 'fdfmethod.lagrange';

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
