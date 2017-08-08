function privvals = absfilterprivvals(this)
%ABSFILTERPRIVVALS   

%   Author(s): R. Losada
%   Copyright 2003 The MathWorks, Inc.

privvals.parent = basefilterprivvals(this);

% Get names of private properties we want to copy
pnames = absfilterprivnames(this);

privvals.this = get(this,pnames);

% [EOF]
