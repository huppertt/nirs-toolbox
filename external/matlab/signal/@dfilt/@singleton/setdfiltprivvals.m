function setdfiltprivvals(this, privvals)
%SETDFILTPRIVVALS   Set private values.

%   Author(s): R. Losada
%   Copyright 2003-2005 The MathWorks, Inc.

setabsfilterprivvals(this, privvals.parent);

% Get names of private properties we want to copy
pnames = dfiltprivnames(this);

set(this,pnames,privvals.this);

% [EOF]
