function setabsfilterprivvals(this, privvals)
%SETABSFILTERPRIVVALS   

%   Author(s): R. Losada
%   Copyright 2003 The MathWorks, Inc.

setbasefilterprivvals(this, privvals.parent);

% Get names of private properties we want to copy
pnames = absfilterprivnames(this);

set(this,pnames,privvals.this);

% [EOF]
