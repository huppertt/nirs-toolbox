function [y,zf] = secfilter(this,x,zi)
%SECFILTER Filter this section.
%   [Y,Zf] = SECFILTER(Hd,X,ZI) filters this section.
%
%   See also DFILT. 

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

Hd = dispatch(this);

reset(Hd); % Reset stages of contained cascade
Hd.PersistentMemory = true;

% Assign States from cascade object to each stage in contained object
count = 0;
for n = 1:nstages(Hd),
    ncoeffs_stage = length(Hd.Stage(n).AllpassCoefficients);
    Hd.Stage(n).States = zi(count+1:count+ncoeffs_stage,:);
    count = count + ncoeffs_stage;
end

if size(x,1) == 1,
    y = filter(Hd,x,1);
else
    y = filter(Hd,x);
end

zf = zeros(size(zi));

% Assign end states back in cascade object
count = 0;
for n = 1:nstages(Hd),
    ncoeffs_stage = length(Hd.Stage(n).AllpassCoefficients);
    zf(count+1:count+ncoeffs_stage,:) = Hd.Stage(n).States;
    count = count + ncoeffs_stage;
end



% [EOF]
