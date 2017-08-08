function nsecs = nsections(Hb)
%NSECTIONS Number of sections in a discrete filter.
%   NSECTIONS(Hb) returns the number of sections in a discrete
%   filter.
%
%   See also DFILT.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

nsecs = base_num(Hb, 'thisnsections');

% [EOF]
