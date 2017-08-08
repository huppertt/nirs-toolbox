function designobj = getdesignobj(~, str, sigonlyflag) %#ok<INUSD>
%GETDESIGNOBJ   Get the designobj.

%   Copyright 1988-2013 The MathWorks, Inc.

if nargin < 2
    str = [];  
end

%#function fmethod.butterbs
designobj.butter = 'fmethod.butterbs';

%#function fmethod.windowbs
designobj.window = 'fmethod.windowbs';

if ~isempty(str)
    designobj = designobj.(str);
end

% [EOF]
