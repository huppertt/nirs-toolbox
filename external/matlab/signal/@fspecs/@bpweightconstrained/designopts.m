function s = designopts(this, dmethod)
%DESIGNOPTS Display the design options.

%   Copyright 2011 The MathWorks, Inc.

hmethod = feval(getdesignobj(this, dmethod));

% Need to pass hspecs so that the design object can know if a band is
% constrained or not and add or remove the weight design options for that
% band accordingly.
s = designopts(hmethod, this);

% [EOF]
