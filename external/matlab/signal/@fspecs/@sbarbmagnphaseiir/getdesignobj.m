function designobj = getdesignobj(~, str)
%GETDESIGNOBJ   Get the design object.

%   Copyright 2005-2010 The MathWorks, Inc.

%#function fmethod.invfreqz2
designobj.iirls = 'fmethod.invfreqz2';

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
