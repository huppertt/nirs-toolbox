function setbasefilterprivvals(this, privvals)
%SETBASEFILTERPRIVVALS   

%   Author(s): R. Losada
%   Copyright 2003 The MathWorks, Inc.

% Get names of private properties we want to copy
pnames = basefilterprivnames(this);

set(this,pnames,privvals);

% [EOF]
