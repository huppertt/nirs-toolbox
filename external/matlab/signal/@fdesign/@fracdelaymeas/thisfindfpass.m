function F = thisfindfpass(this,hfilter,cellarray)
%THISFINDFPASS  If both Fpass and fderr are empty we cannot find an Fpass. 

%   Author(s): V. Pellissier
%   Copyright 2006 The MathWorks, Inc.

% If both Fpass and fderr are empty, use FracDelayError to find Fpass.

N = 2^10;
[gpd w fderr updwn Fs] = deal(cellarray{:});
updwn = strcmpi(updwn,'up');

if isempty(fderr),
    F = [];
else
    nomgdp = floor(order(hfilter)/2);
    fd = hfilter.FracDelay;
    idealgpd = nomgdp + fd;
    if ~this.NormalizedFrequency,
        % always work in samples for the group delay
        fderr = fderr*Fs; 
    end
    [w_lo, w_hi] = lclfindfpass(w, gpd-idealgpd, fderr, updwn,[0 Fs/2]); 
    if isempty(w_lo) || isempty(w_hi)
        F = [];
        return;
    end
    % Refine search
    [gpd, w] = grpdelay(hfilter, linspace(w_lo, w_hi, N), Fs);
    [w_lo, w_hi] = lclfindfpass(w, gpd-idealgpd, fderr, updwn);
    if isempty(w_lo) || isempty(w_hi)
        F = [];
        return;
    end
    F = (w_lo+w_hi)/2;
end

% -------------------------------------------------------------------------
function [w_lo, w_hi] = lclfindfpass(w, gpd, fderr, up, Frange)

% Remove the areas to "ignore".
if nargin > 4
    lo_indx = find(w < Frange(1));
    hi_indx = find(w > Frange(2));
    w([lo_indx hi_indx]) = [];
    gpd([lo_indx hi_indx]) = [];
end

if up
    indx = find(gpd > -fderr, 1, 'first');
    w_lo = w(max(1, indx-1));
    w_hi = w(indx);
else
    indx = find(gpd > -fderr, 1, 'last');
    w_lo = w(indx);
    w_hi = w(min(indx+1, length(w)));
end

% [EOF]
