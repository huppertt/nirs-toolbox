function designobj = getdesignobj(this, str)
%GETDESIGNOBJ Get the design object
%   OUT = GETDESIGNOBJ(ARGS) <long description>

%   Copyright 2008 The MathWorks, Inc.

%#function fmethod.sqrtrcoswin
designobj.window = 'fmethod.sqrtrcoswin';

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
