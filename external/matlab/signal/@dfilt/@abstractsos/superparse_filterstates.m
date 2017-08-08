function hTar = superparse_filterstates(Hd,hTar)
%SUPERPARSE_FILTERSTATES Store filter states in hTar for df1sos and df1tsos

%   Copyright 2009 The MathWorks, Inc.

% Extract current filter states
IC = getinitialconditions(Hd);

% If the MapStates stored in hTar is not 'on', set the initial condition to
% 0.
if ~strcmpi(hTar.MapStates,'on')
    IC.Num = zeros(size(IC.Num));
    IC.Den = zeros(size(IC.Den));
end

% Store the filter states
setprivstates(hTar,IC);

% [EOF]
