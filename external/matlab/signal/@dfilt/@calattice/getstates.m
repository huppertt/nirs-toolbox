function S = getstates(Hm,S)
%GETSTATES Overloaded get for the States property.

% This should be a private method

%   Author: R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

S = Hm.HiddenStates;
if nstates(Hm)==0,
     S = [];
     return
end
n1 = 0;
ncoeffs = Hm.ncoeffs;
if ~isempty(ncoeffs), n1 = Hm.ncoeffs(1); end
S = [S(1:n1,:);S(n1+2:end-1,:)];
