function S = thissetstates(Hd,S)
%THISSETSTATES Overloaded set for the States property.

% This should be a private method

%   Author: V. Pellissier, P. Costa
%   Copyright 1988-2004 The MathWorks, Inc.

% Hd.ncoeffs is being set to [0 1] even though the constructor sets it
% to [1 1].
if isempty(Hd.ncoeffs)
    return;
end

if  Hd.ncoeffs(1)== 0, Hd.ncoeffs(1) = 1; end
nb = Hd.ncoeffs(1);
if length(Hd.ncoeffs) ==1, Hd.ncoeffs(2) = 1; end
na = Hd.ncoeffs(2);
Snum = [];
Sden = [];

% Check data type, quantize if needed
S = validatestatesobj(Hd.filterquantizer, S);

% If the states are not empty, convert linear to circular states
if ~isempty(S),

    % Convert linear to circular states
    tapIndex = Hd.tapindex;

    % In R13, tapIndex was a scalar.
    if length(tapIndex == 1), tapIndex = [tapIndex 0]; end

    % NOTE: The DF1FILTER.DLL filter implementation expects the number of
    % numerator states to be NB whereas the number of denominator states
    % should be NA-1, due to the use of circular buffers.

    % Because the spec for setting the states for the DF1, DF1T, DF1SOS, and
    % DF1TSOS allows one to set both a FILTSTATES.DFIIR and a double
    % vector, we create a FILTSTATES.DFIIR object in ziscalarexpand (which is
    % always called before thissetstates -- see @abstractfilter/schema.m)

    % Create numerator and denominator states with the correct
    % attributes; Circular buffer expects a flipped version of the
    % states.
    Snum_lin = flipud(S.Numerator);

    % Add an extra zero to the numerator states
    Snum_lin = prependzero(Hd.filterquantizer,Snum_lin);
    Sden_lin = flipud(S.Denominator);
    Snum = S.Numerator;
    Sden = S.Denominator;

    %
    % Numerator states
    %
    zNIdx = tapIndex(1); zNIdx = double(zNIdx);

    % Extract the numerator portion of the states
    Snum(zNIdx+1:nb,:) = Snum_lin(1:nb-zNIdx,:);
    Snum(1:zNIdx,:) = Snum_lin(nb-zNIdx+1:nb,:);

    %
    % Denominator States
    %
    zDIdx = tapIndex(2); zDIdx = double(zDIdx);
    Sden(zDIdx+1:end,:) = Sden_lin(1:na-1-zDIdx,:);
    Sden(1:zDIdx,:) = Sden_lin(na-zDIdx:end,:);

    % Update the FILTSTATES object.
    S.Numerator = Snum;
    S.Denominator = Sden;
else
    % The DLL expects at least one state (of zero) for the numerator side,
    % therefore, we create a null state for the numerator and create a
    % corresponding for the denominator so that ZIEXPAND doesn't
    % error.
    S.Numerator   = nullnumstate(Hd.filterquantizer);
    S.Denominator = nulldenstate(Hd.filterquantizer);
end

% The DFIIRSTATES object contains the circular states (ready for the
% filter implementation)
Hd.HiddenStates = S;
