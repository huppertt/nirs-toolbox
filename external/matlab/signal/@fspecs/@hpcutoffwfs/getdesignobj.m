function designobj = getdesignobj(this, str)
%GETDESIGNOBJ   Get the design object.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

%#function fmethod.cheby2hp
designobj.cheby2 = 'fmethod.cheby2hp';

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
