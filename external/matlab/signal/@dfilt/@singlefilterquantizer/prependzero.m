function S = prependzero(q,S)
%PREPENDZERO   

%   Author(s): V. Pellissier
%   Copyright 1999-2003 The MathWorks, Inc.

S = [single(zeros(1,size(S,2)));S];


% [EOF]
