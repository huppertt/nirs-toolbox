function designobj = getdesignobj(~, str, sigonlyflag) %#ok<INUSD>
%GETDESIGNOBJ   Get the designobj.

%   Copyright 1988-2013 The MathWorks, Inc.

if nargin < 2
    str = [];  
end

%#function fmethod.cheby2hp
designobj.cheby2 = 'fmethod.cheby2hp';

if ~isempty(str)
    designobj = designobj.(str);
end

% [EOF]
