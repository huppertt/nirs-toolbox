function [F, A] = getmask(this)
%GETMASK Get the mask.

%   Copyright 2005-2011 The MathWorks, Inc.

F = this.Frequencies;
A = this.FreqResponse;

if ~this.NormalizedFrequency
  F = F/(this.Fs/2);
end

% [EOF]
