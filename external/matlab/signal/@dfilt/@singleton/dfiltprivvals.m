function privvals = dfiltprivvals(this)
%DFILTPRIVVALS   

%   Author(s): R. Losada
%   Copyright 2003-2005 The MathWorks, Inc.

privvals.parent = absfilterprivvals(this);

% Get names of private properties we want to copy
pnames = dfiltprivnames(this);

privvals.this = get(this,pnames);


% [EOF]
