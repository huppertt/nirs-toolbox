function c = cscalefactors(h,opts)
%CSCALEFACTORS   Cumulative scale factors

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

nb = nsections(h);

c = zeros(1,nb+1);

c = lclcscalefactors(h,c,nb,opts);

% [EOF]
