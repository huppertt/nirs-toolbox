function n = nadd(this)
%NADD Returns the number of adders  

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

refsosMatrix = this.refsosMatrix;
n = length(find(refsosMatrix~=0))-2*nsections(this);

% [EOF]
