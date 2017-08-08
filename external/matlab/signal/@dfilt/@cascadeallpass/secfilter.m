function [y,zf] = secfilter(this,x,zi)
%SECFILTER Filter this section.
%   [Y,Zf] = SECFILTER(Hd,X,ZI) filters this section.
%
%   See also DFILT. 

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

Hd = dispatch(this);

% Determine index mappings for states of all stages
nst = nstages(Hd);

% Special case the first stage
ncoeffs{1} = length(Hd.Stage(1).AllpassCoefficients);

s(1).ceil  = 1:ncoeffs{1};
s(1).floor = ncoeffs{1}+1:2*ncoeffs{1};

% Now find indexes for remaining stages
for k = 2:nst,
    ncoeffs{k} = length(Hd.Stage(k).AllpassCoefficients);
    s(k).ceil  = s(k-1).floor(1):s(k-1).floor(1)+ncoeffs{k}-1;
    % The floor index depends on whether or not the previous stage had more
    % coeffs
    floorstartidx = max(s(k-1).floor(end),s(k).ceil(end))+1;
    s(k).floor = floorstartidx:floorstartidx+ncoeffs{k}-1;
end

reset(Hd); % Reset stages of contained cascade
Hd.PersistentMemory = true;

% Map states using indexes
for k = 1:nst,
    Hd.Stage(k).States = [zi(s(k).ceil,:);zi(s(k).floor,:)];
end


if size(x,1) == 1,
    y = filter(Hd,x,1);
else
    y = filter(Hd,x);
end


zf = zeros(size(zi));

% Assign end states back in cascade object
for k = 1:nst,
    zf(s(k).ceil,:)  = Hd.Stage(k).States(1:ncoeffs{k},:);
    zf(s(k).floor,:) = Hd.Stage(k).States(ncoeffs{k}+1:end,:);
end

% [EOF]
