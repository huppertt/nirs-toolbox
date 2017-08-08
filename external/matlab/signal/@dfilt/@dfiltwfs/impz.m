function [y, t] = impz(Hd, varargin)
%IMPZ Returns the impulse response

%   Author(s): J. Schickler
%   Copyright 1988-2009 The MathWorks, Inc.

[y, t] = timeresp(Hd, @lclimpz, varargin{:});

% -----------------------------------------------------------
function [y, t] = lclimpz(G, N, Fs)

if isempty(Fs), Fs = 1; end
[y, t] = impz(G, N, Fs);

% [EOF]
