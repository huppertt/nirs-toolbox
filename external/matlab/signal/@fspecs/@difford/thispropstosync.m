function p = thispropstosync(this, p)
%THISPROPSTOSYNC   

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

% Remove FilterOrder because it must be even for certain specs and odd for
% others.
p(strmatch('FilterOrder',p)) = [];

% [EOF]
