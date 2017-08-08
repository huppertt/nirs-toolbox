function F = findfpass(this, hfilter, Fpass, Apass, direction, Frange, idealfcn)
%FINDFPASS   Find the FPass value.
%   F = FINDFPASS(H, HD, FP, AP, DIR, PLACE) returns the Fpass value in F
%   for the given filter object HD.  FP is the "known" Fpass (if it
%   exists), AP is the "known" Apass, DIR is the string 'up' or 'down'
%   which indicates whether the response is moving up from left to right or
%   down.  PLACE is the string 'first' or 'last' which indicates if this is
%   the first passband or the last.  If FP is not [], the algorithm will
%   simply return it.  If AP is [] Fpass cannot be determined.

%   Author(s): J. Schickler
%   Copyright 2005-2009 The MathWorks, Inc.

if ~isempty(Fpass)
    
    % If we are passed Fpass simply return it, do not measure.
    F = Fpass;
elseif ~isempty(Apass)
    
    % If we are passed Apass, but not Fpass, we can measure.
    
    % Set up the default NFFT.
    N = 2^10;

    if strcmpi(direction, 'up'), up = true;
    else,                        up = false; end

    if this.NormalizedFrequency, Fs = 2;
    else,                        Fs = this.Fs; end
        
    % Perform a FREQZ which we will measure.
    [h, w] = freqz(hfilter, N, Fs);
    
    h = abs(h);
    
    if nargin < 6
        Frange = [0 Fs/2];
    end
    
    if nargin > 6
        if iscell(idealfcn)
            idealh = feval(idealfcn{1}, w/(Fs/2), idealfcn{2:end});
        else
            idealh = feval(idealfcn, w/(Fs/2))
        end
        h = h-idealh+1;
        h(h < 0) = 0;
    end

    % Convert the FREQZ to db magnitude and normalize the highest value.
    h = db(h);

    hmax = max(h);
    h = h-hmax;

    [w_lo, w_hi] = lclfindfpass(w, h, Apass, up, Frange);
    
    if isempty(w_lo) || isempty(w_hi)
        F = [];
        return;
    end

    [h, w] = freqz(hfilter, linspace(w_lo, w_hi, N), Fs);
    h      = abs(h);

    if nargin > 6
        if iscell(idealfcn)
            idealh = feval(idealfcn{1}, w/(Fs/2), idealfcn{2:end});
        else
            idealh = feval(idealfcn, w/(Fs/2))
        end
        h = h-idealh+1;
    end
    
    % Convert to dB and use the previously determined hmax to normalize
    % because the new 'h' vector does not include it.
    h = db(h);
    h = h-hmax;

    [w_lo, w_hi] = lclfindfpass(w, h, Apass, up, Frange);
    
    if isempty(w_lo) || isempty(w_hi)
        F = [];
        return;
    end

    F = (w_lo+w_hi)/2;

else
    if nargin<7,
        idealfcn = [];
    end
    F = thisfindfpass(this,hfilter,idealfcn);
end

% -------------------------------------------------------------------------
function [w_lo, w_hi] = lclfindfpass(w, h, Apass, up, Frange)

% Remove the areas to "ignore".
if nargin > 4
    lo_indx = find(w < Frange(1));
    hi_indx = find(w > Frange(2));
    w([lo_indx;hi_indx]) = [];
    h([lo_indx;hi_indx]) = [];
end

if up
    indx = find(h > -Apass, 1, 'first');
    w_lo = w(max(1, indx-1));
    w_hi = w(indx);
else
    indx = find(h > -Apass, 1, 'last');
    w_lo = w(indx);
    w_hi = w(min(indx+1, length(w)));
end

% [EOF]
