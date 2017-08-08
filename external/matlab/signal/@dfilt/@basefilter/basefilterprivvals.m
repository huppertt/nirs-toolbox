function privvals = basefilterprivvals(this)
%BASEFILTERPRIVVALS   

%   Author(s): R. Losada
%   Copyright 2003 The MathWorks, Inc.

% Get names of private properties we want to copy
pnames = basefilterprivnames(this);

privvals = get(this,pnames);

% [EOF]
