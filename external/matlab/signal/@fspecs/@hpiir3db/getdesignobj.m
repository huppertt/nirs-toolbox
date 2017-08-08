function designobj = getdesignobj(~, str, sigonlyflag) %#ok<INUSD>
%GETDESIGNOBJ Get the design object.

%   Copyright 2011-2013 The MathWorks, Inc.

if nargin < 2
    str = [];  
end

%#function fmethod.iirmaxflathp
designobj.butter = 'fmethod.iirmaxflathp';

if ~isempty(str)
    designobj = designobj.(str);
end

% [EOF]
