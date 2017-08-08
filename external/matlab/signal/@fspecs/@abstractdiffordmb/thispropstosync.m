function p = thispropstosync(~, p)
%THISPROPSTOSYNC   

%   Copyright 2005-2011 The MathWorks, Inc.

% Remove FilterOrder because it must be even for certain specs and odd for
% others.
p(strcmp('FilterOrder',p)) =[];

% [EOF]
