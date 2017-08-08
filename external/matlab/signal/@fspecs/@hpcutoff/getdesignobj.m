function designobj = getdesignobj(~, str, sigonlyflag) %#ok<INUSD>
%GETDESIGNOBJ   Get the designobj.

%   Copyright 1988-2013 The MathWorks, Inc.

if nargin < 2
    str = [];  
end

%#function fmethod.butterhp
designobj.butter = 'fmethod.butterhp';

%#function fmethod.windowhp
designobj.window = 'fmethod.windowhp';

if ~isempty(str)
    designobj = designobj.(str);
end

% [EOF]
