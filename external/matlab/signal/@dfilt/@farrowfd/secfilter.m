function [y,z] = secfilter(this,x,d,z)
%SECFILTER   

%   Author(s): V. Pellissier
%   Copyright 2005-2006 The MathWorks, Inc.
C = this.privcoeffs;
C = C.';
[y,z] = farrowfdfilter(this.filterquantizer,C,x,d,z);

% [EOF]
