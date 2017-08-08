function ic = getinitialconditions(Hd)
%GETINITIALCONDITIONS Get the initial conditions

%   Copyright 2009 The MathWorks, Inc.

s    = double(Hd.States);
nsts  = size(Hd.sosMatrix,1)*2;
nchan = prod(size(s))/nsts;
ic    = reshape(s,nsts,nchan);

% [EOF]
