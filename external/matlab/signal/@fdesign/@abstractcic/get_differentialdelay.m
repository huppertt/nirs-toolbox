function diffd = get_differentialdelay(this, diffd)
%GET_DIFFERENTIALDELAY   PreGet function for the 'differentialdelay'
%property.

%   Author(s): P. Costa
%   Copyright 2005 The MathWorks, Inc.

diffd = this.privDifferentialDelay;
if isempty(diffd), diffd = 1; end


% [EOF]
