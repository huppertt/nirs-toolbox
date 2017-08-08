function S = getstates(Hd,S)
%GETSTATES Overloaded get for the States property.

% This should be a private method

%   Author: V. Pellissier
%   Copyright 1988-2004 The MathWorks, Inc.

% Convert circular states to linear states
Scir = Hd.HiddenStates;

if isempty(Scir)
    S = [];
else

    tapIndex = Hd.tapIndex(1)+1; % Due to the 1-based indexing
    S = [Scir(tapIndex+1:end,:); Scir(1:tapIndex-1,:)];
end

% [EOF]
