function [F, A] = getmask(this)
%GETMASK Get the mask.

%   Copyright 2005-2011 The MathWorks, Inc.

p = propstoadd(this);
% Remove NormalizedFrequency, Fs, FilterOrder and Nbands
p([1 2 3 4]) = [];

F = [];
A = [];
for i=1:2:length(p),
    F = [F this.(p{i})]; %#ok<*AGROW>
    A = [A this.(p{i+1})];
end

if ~this.NormalizedFrequency
  F = F/(this.Fs/2);
end

% [EOF]
