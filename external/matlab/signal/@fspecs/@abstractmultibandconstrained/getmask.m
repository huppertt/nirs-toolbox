function [F, A] = getmask(this)
%GETMASK Get the mask.

%   Copyright 2011 The MathWorks, Inc.

F = [];
A = [];
for i=1:this.NBands
    F = [F this.(sprintf('%s%d%s','B',i,'Frequencies'))]; %#ok<*AGROW>
    A = [A this.(sprintf('%s%d%s','B',i,'Amplitudes'))];
end

if ~this.NormalizedFrequency
  F = F/(this.Fs/2);
end

% [EOF]
