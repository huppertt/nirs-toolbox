function b = isvalidparallel(this)
%ISVALIDPARALLEL   True if the object is validparallel.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

b = true;
if nstages(this) > 0,
    firstrcf = prod(getratechangefactors(this.stage(1)),1);
    firstrcf = firstrcf(1)/firstrcf(2);
    count = 2;
    while b && count <= nstages(this),
        stagercf = prod(getratechangefactors(this.stage(count)),1);
        stagercf = stagercf(1)/stagercf(2);
        b = b && all(firstrcf == stagercf);
        count = count + 1;
    end
end


% [EOF]
