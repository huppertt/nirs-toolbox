function F = findfrequency(this, hfilter, A, direction, place)
%FINDFREQUENCY   Find the frequency point for the given amplitude
%   F = FINDFREQUENCY(H, HD, A, DIR, PLACE) returns the frequency point F
%   for the amplitude A in the filter object HD.  A is in linear units and
%   is relative to the nominal gain of the filter.
%
%   F = FINDFREQUENCY(H, 


%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

error(nargchk(4,5,nargin,'struct'));

N = 2^10;

ngain = nominalgain(hfilter);
if ~isempty(ngain)
    A = ngain*A;
end

if this.NormalizedFrequency, Fs = 2;
else,                        Fs = this.Fs; end

if strcmpi(direction, 'up'), up = true;
else,                        up = false; end

if strcmpi(place, 'first'), first = true;
else,                       first = false; end

% find first point that fails a3db
[h,w] = freqz(hfilter, N, Fs);
h = abs(h);

[w_lo, w_hi] = lclfindfrequency(w, h, A, up, first);

[F htest] = getF(hfilter, [w_lo w_hi], Fs);

count = 1;

while abs(A-htest) > 1e-5 && count < 4

    % Do a second refinement if we did not hit the point exactly.

    [h, w] = freqz(hfilter, linspace(w_lo, w_hi, N), Fs);
    h = abs(h);

    [w_lo, w_hi] = lclfindfrequency(w, h, A, up, first);

    [F htest] = getF(hfilter, [w_lo w_hi], Fs);
    
    count = count+1;
end

% -------------------------------------------------------------------------
function [w_lo, w_hi] = lclfindfrequency(w, h, A, up, first)

if first
    if up
        indx = find(h >= A, 1, 'first');
    else
        indx = find(h <= A, 1, 'first');
    end
    if isempty(indx), indx = length(w); end

    w_lo = w(max(1, indx-1));
    w_hi = w(indx);
else
    if up
        indx = find(h <= A, 1, 'last');
    else
        indx = find(h >= A, 1, 'last');
    end
    if isempty(indx), indx = 1; end
    w_lo = w(indx);
    w_hi = w(min(length(w), indx+1));
end


% -------------------------------------------------------------------------
function [F, htest] = getF(hfilter, w, Fs)

F = mean(w);

htest = freqz(hfilter, [0 F], Fs);
htest = abs(htest(2));

% [EOF]
