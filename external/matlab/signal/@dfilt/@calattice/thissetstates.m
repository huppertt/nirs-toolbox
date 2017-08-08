function S = thissetstates(Hm,S)
%THISSETSTATES Overloaded set for the States property.

% This should be a private method

%   Author: R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

o = zeros(1,size(S,2));
n1 = 0;
ncoeffs = Hm.ncoeffs;
if ~isempty(ncoeffs), n1 = Hm.ncoeffs(1); end
S = [S(1:n1,:);o;S(n1+1:end,:);o];
if isempty(S), S=[0;0]; end

Hm.HiddenStates = S;

