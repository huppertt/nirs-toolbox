function [Hm,meas] = multisection(this,M,R)
%CICDESIGN   Shared design gateway for CICDECIM/CICINTERP.

%   Author(s): P. Costa
%   Copyright 2005 The MathWorks, Inc.

do = this.getdesignobj('multisection');

Hdo = feval(do);

Hm = design(Hdo,this,R,M);

% [EOF]
