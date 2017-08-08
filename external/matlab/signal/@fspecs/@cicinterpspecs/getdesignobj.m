function designobj = getdesignobj(this, str)
%GETDESIGNOBJ   Get the designobj.

%   Author(s): P. Costa
%   Copyright 2005 The MathWorks, Inc.

%#function fdfmethod.designcicinterp
designobj.multisection = 'fdfmethod.designcicinterp';

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
