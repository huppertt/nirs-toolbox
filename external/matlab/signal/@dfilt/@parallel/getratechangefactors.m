function rcf = getratechangefactors(this)
%GETRATECHANGEFACTORS   Get the ratechangefactors.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

checkvalidparallel(this);

if nstages(this) > 0,
    rcf = prod(getratechangefactors(this.stage(1)),1);
else
    rcf = [1 1];
end


% [EOF]
