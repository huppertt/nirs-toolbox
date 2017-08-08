function hTar = parse_filterstates(Hd,hTar)
%PARSE_FILTERSTATES Store filter states in hTar for realizemdl

%   Copyright 2009 The MathWorks, Inc.

% Extract current filter states
IC = getinitialconditions(Hd);

% If the MapStates stored in hTar is not 'on', set the initial condition to
% 0.
if ~strcmpi(hTar.MapStates,'on')
    IC = zeros(size(IC));
end

% Store the filter states
setprivstates(hTar,IC);


% [EOF]
