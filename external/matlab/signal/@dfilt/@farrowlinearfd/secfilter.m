function [y,z] = secfilter(this,x,d,z)
%SECFILTER   

%   Author(s): V. Pellissier
%   Copyright 2005-2006 The MathWorks, Inc.

q = this.filterquantizer;
[y,z] = linearfdfilter(q,x,d,z);


% [EOF]
