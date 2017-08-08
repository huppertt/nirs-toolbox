function minfo = measureinfo(this)
%MEASUREINFO   

%   Copyright 2010 The MathWorks, Inc.

p = propstoadd(this);

minfo.NomGrpDelay = this.NomGrpDelay; % this value is computed at design time

% Remove NormalizedFrequency, Fs, FilterOrder and Nbands
p([1 2 3 4]) = [];

F = [];
Gd = [];
for i=1:2:length(p),
    F = [F this.(p{i})]; %#ok<*AGROW>
    Gd = [Gd this.(p{i+1})];
end

minfo.Frequencies = F;
minfo.GroupDelay = Gd; 

% [EOF]
