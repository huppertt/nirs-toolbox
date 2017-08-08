function designobj = getdesignobj(~, str, sigonlyflag) %#ok<INUSD>
%GETDESIGNOBJ   Get the designobj.

%   Copyright 1988-2013 The MathWorks, Inc.

if nargin < 2
    str = [];  
end

%#function fmethod.ellipbpastop
designobj.ellip = 'fmethod.ellipbpastop';

if ~isempty(str)
    designobj = designobj.(str);
end

% [EOF]
