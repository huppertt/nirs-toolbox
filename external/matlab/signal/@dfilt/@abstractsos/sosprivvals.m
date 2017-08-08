function privvals = sosprivvals(this)
%SOSPRIVVALS   

%   Author(s): R. Losada
%   Copyright 2003 The MathWorks, Inc.

privvals.parent = dfiltprivvals(this);

% Get names of private properties we want to copy
pnames = sosprivnames(this);

privvals.this = get(this,pnames);


% [EOF]
