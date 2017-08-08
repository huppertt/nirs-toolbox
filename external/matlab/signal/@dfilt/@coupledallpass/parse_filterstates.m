function hTar = parse_filterstates(Hd,hTar)
%PARSE_FILTERSTATES Store filter states in hTar for realizemdl

%   Copyright 2009 The MathWorks, Inc.

% Extract current filter states
ic = getinitialconditions(Hd);

% If the MapStates stored in hTar is not 'on', set the initial condition to
% 0.
if ~strcmpi(hTar.MapStates,'on')
    ic = zeros(size(ic));
end

% Separate states for each allpass stage
nstates1 = length(Hd.Allpass1);
IC{1} = ic(1:nstates1,:);
IC{2} = ic(nstates1+1:end,:);

% Store the filter states
setprivstates(hTar,IC);
% [EOF]
