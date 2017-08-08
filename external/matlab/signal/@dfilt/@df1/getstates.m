function S = getstates(Hd,S)
%GETSTATES Overloaded get for the States property.

% This should be a private method

%   Author: V. Pellissier
%   Copyright 1988-2004 The MathWorks, Inc.

ncoeffs = get(Hd, 'ncoeffs');
if isempty(ncoeffs)
    return;
end

nb = ncoeffs(1);
na = ncoeffs(2);

% HiddenStates contains the states object which has the circular
% representation of both Numerator and Denominator states.
S = copy(Hd.HiddenStates);

% Convert the circular to linear states
if ~isempty(S),

    % Extract the circular states for both the numerator and denominator
    Snum_cs = S.Numerator;
    Sden_cs = S.Denominator;

    % Create default vectors (Snum & Sden) of the right class to contain
    % the linear states.  We are making copies of the *_cs vectors here to
    % preserve the data type of the states.
    Snum = Snum_cs([]);
    Sden = Sden_cs([]);

    tapIndex = Hd.tapIndex;
    % In R13, tapIndex was a scalar.
    if length(tapIndex == 1), tapIndex = [tapIndex 0]; end

    % Numerator states
    zfIdx = tapIndex(1);
    if zfIdx<1
        zfIdx = nb;
    end
    for j = 1:(nb-1)
        Snum(j,:) = Snum_cs(zfIdx,:);
        zfIdx = zfIdx-1;
        if zfIdx<1
            zfIdx = nb;
        end
    end

    % Denominator States
    zfIdx = tapIndex(2);
    if zfIdx<1
        zfIdx = na-1;
    end
    for j = 1:na-1
        Sden(j,:) = Sden_cs(zfIdx,:);
        zfIdx = zfIdx-1;
        if zfIdx<1
            zfIdx = na-1;
        end
    end

    % Update the properties of the states object (with the linear states)
    S.Numerator = Snum;
    S.denominator = Sden;
end
