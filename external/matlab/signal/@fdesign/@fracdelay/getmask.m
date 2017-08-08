function [F, A, Gd] = getmask(this, fcns, rcf, specs)
%GETMASK   Get the mask.

%   Author(s): V. Pellissier
%   Copyright 2005-2006 The MathWorks, Inc.

F = [0 1]*fcns.getfs()/2;
A = fcns.formatapass(0);
Gd = floor(this.FilterOrder/2) + getfracdelay(this.CurrentSpecs);
Gd = [Gd Gd];

% [EOF]
