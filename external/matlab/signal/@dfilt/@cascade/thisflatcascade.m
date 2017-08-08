function Hflat = thisflatcascade(this,Hflat)
%THISFLATCASCADE Add singletons to the flat list of filters Hflat 

%   Copyright 2008 The MathWorks, Inc.

N = nstages(this);
for i=1:N,
    Hflat = thisflatcascade(this.Stage(i),Hflat);
end
