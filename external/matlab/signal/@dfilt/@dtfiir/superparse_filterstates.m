function hTar = superparse_filterstates(Hd,hTar)
%SUPERPARSE_FILTERSTATES Store filter states in hTar for df1sos and df1tsos

%   Copyright 2009 The MathWorks, Inc.

% Extract current filter states
IC = getinitialconditions(Hd);

% Make states the same order. When OptimizeZeros is 'off', realizemdl
% performs zero-padding so that the numerator and denominator are balance.
% We, thus, need to zero-pad the filter states as well.
max_order = max(length(IC.Num),length(IC.Den));
IC.Num = makemaxorder(IC.Num,max_order);
IC.Den = makemaxorder(IC.Den,max_order);

% If the MapStates stored in hTar is not 'on', set the initial condition to
% 0.
if ~strcmpi(hTar.MapStates,'on')
    IC.Num = zeros(size(IC.Num));
    IC.Den = zeros(size(IC.Den));
end

% Store the filter states
setprivstates(hTar,IC);

%--------------------------------------------------------------------------
function ic = makemaxorder(ic,maxorder)

currentorder = length(ic);
M = abs(maxorder-currentorder);
if currentorder < maxorder
    % If currentorder is less than maxorder, the zero padding is required
    % so that the filter states are balance.
    ic = [ic; zeros(M,1)];
elseif currentorder > maxorder
    % If the currentorder (coefficient length) is larger than maxorder,
    % then remove the exceeding states as they are zeros.
    ic = ic(1:end-M);
end


% [EOF]
