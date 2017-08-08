function s = designopts(this, dmethod, sigonlyflag)
%DESIGNOPTS   Display the design options.

%   Copyright 2005-2013 The MathWorks, Inc.

if nargin == 3
    hmethod = feval(getdesignobj(this, dmethod, sigonlyflag));
else
    hmethod = feval(getdesignobj(this, dmethod));
end

s = designopts(hmethod);

% [EOF]
