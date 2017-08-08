function designobj = getdesignobj(~, str, sigonlyflag) %#ok<INUSD>
%GETDESIGNOBJ   Get the design object.

%   Copyright 2005-2013 The MathWorks, Inc.

if nargin < 2
    str = [];  
end

%#function fmethod.butterbs
designobj.butter = 'fmethod.butterbs';

if ~isempty(str)
    designobj = designobj.(str);
end

% [EOF]
