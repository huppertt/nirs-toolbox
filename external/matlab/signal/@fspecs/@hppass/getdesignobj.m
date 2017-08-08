function designobj = getdesignobj(~, str, sigonlyflag) %#ok<INUSD>
%GETDESIGNOBJ   Get the designobj.

%   Copyright 1988-2013 The MathWorks, Inc.

if nargin < 2
    str = [];  
end

%#function fmethod.cheby1hp
designobj.cheby1 = 'fmethod.cheby1hp';

if ~isempty(str)
    designobj = designobj.(str);
end

% [EOF]
