function designobj = getdesignobj(~, str, sigonlyflag) %#ok<INUSD>
%GETDESIGNOBJ   Get the design object.

%   Copyright 2008-2013 The MathWorks, Inc.

if nargin < 2
    str = [];  
end

%#function fmethod.firclsbp
designobj.fircls = 'fmethod.firclsbp';

if ~isempty(str)
    designobj = designobj.(str);
end

% [EOF]
