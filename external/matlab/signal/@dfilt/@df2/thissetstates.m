function S = thissetstates(Hd,S)
%THISSETSTATES Overloaded set for the States property.

% This should be a private method

%   Author: V. Pellissier
%   Copyright 1988-2004 The MathWorks, Inc.

%   Zero-pads initial states of the discrete-time filter Hd before
%   passing to the C++ filtering function.  Zero padding of the initial
%   states is done to prevent special-case code for the first and last
%   states.

% The circular buffer implementation in DF2FILTER.DLL expects the number of
% states to be MAX(Order_Num,Order_Num) + 1;

if isempty(S),
    S = nullstate1(Hd.filterquantizer);
else
    % Check data type, quantize if needed
    S = validatestates(Hd.filterquantizer, S);
    % Prepended the state vector with one extra zero per channel
    S = prependzero(Hd.filterquantizer, S);
end

% If there is no tap index, we are in an init state, S is already correct.
if isempty(Hd.tapindex)
    return;
end

% Convert linear to circular states
tapIndex = Hd.tapindex(1);

[nr ncol] = size(S);
Scir = S;
Scir(tapIndex+1:end,:) = S(1:nr-tapIndex,:);
Scir(1:tapIndex,:) = S(nr-tapIndex+1:end,:);

Hd.HiddenStates = Scir;

% [EOF]
