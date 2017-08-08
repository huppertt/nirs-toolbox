function [F, A] = getmask(this, fcns, ~, ~)
%GETMASK Get the mask.

%   Copyright 2005-2011 The MathWorks, Inc.

w = warning('off'); %#ok<WNOFF>
[F, A] = getmask(this.CurrentSpecs);
A = fcns.getarbmag(A);
F = F*(fcns.getfs()/2); 
warning(w);

% [EOF]
