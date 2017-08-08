function c = cost(this, archobj)
%COST   

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

c = [];
if nargin < 2,
    for i=1:length(this),
        c = [c evalcost(this(i))];
    end
end

% [EOF]
