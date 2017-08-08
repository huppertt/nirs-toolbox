function F = findfstop(this, hfilter, Fstop, Astop, direction, varargin)
%FINDFSTOP   Find the FStop value.

%   Author(s): J. Schickler
%   Copyright 2005-2009 The MathWorks, Inc.

if ~isempty(Fstop)
    F = Fstop;
elseif ~isempty(Astop)
    N = 2^10;

    if this.NormalizedFrequency, Fs = 2;
    else,                        Fs = this.Fs; end
    
    % Get the nominal gain
    gain = nominalgain(hfilter);
    if isempty(gain)
        gain = 1;
    end

    if strcmpi(direction, 'up'), up = true;
    else,                        up = false; end

    % Perform a FREQZ and normalize by the nominal gain.
    [h, w] = freqz(hfilter, N, Fs);
    h      = db(gain)-db(abs(h))+eps^(1/3);
    
    [w_lo, w_hi] = lclfindfstop(w, h, Astop, up, varargin{:});
    
    if isempty(w_lo) || isempty(w_hi)
        F = [];
    else

        % Perform a 2nd refinement.
        [h, w] = freqz(hfilter, linspace(w_lo, w_hi, N), Fs);
        h      = db(gain)-db(abs(h));

        [w_lo, w_hi] = lclfindfstop(w, h, Astop, up, varargin{:});

        F = mean([w_lo w_hi]);
    end
else
    F = [];
end

% -------------------------------------------------------------------------
function [w_lo, w_hi] = lclfindfstop(w, h, Astop, up, Frange)

% Remove the areas to "ignore".
if nargin > 4
    lo_indx = find(w < Frange(1));
    hi_indx = find(w > Frange(2));
    w([lo_indx; hi_indx]) = [];
    h([lo_indx; hi_indx]) = [];
end

if up
    indx = find(h > Astop, 1, 'last');
    w_lo = w(indx);
    w_hi = w(min(indx+1, length(w)));
else
    indx = find(h > Astop, 1, 'first');
    w_lo = w(max(1, indx-1));
    w_hi = w(indx);
end

% [EOF]
