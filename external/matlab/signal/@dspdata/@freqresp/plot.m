function varargout = plot(this)
%PLOT

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

% Overloaded to do the "dual plot".
hax = newplot;

% Create the Magnitude Plot
hax = gca;

set(hax, 'YAxisLocation', 'Left');

w = get(this, 'Frequencies');
h = get(this, 'Data');

if this.NormalizedFrequency
    xlbl = getfreqlbl('rad');
    w    = w/pi;
else
    [w, m, xunits] = engunits(w);
    xlbl = getfreqlbl([xunits 'Hz']);
end

hl = line(w, 20*log10(abs(h)), 'Parent', hax);

ylabel(hax, getString(message('signal:dspdata:dspdata:MagnitudedB')));
xlabel(hax, xlbl);
title(hax, getString(message('signal:dspdata:dspdata:MagnitudeAndPhaseResponse')));

% Add a second axes in the same location as GCA
hax2 = axes('Units', get(gca, 'Units'), ...
    'Position', get(gca, 'Position'), ...
    'YAxisLocation', 'Right', ...
    'Color', 'none');

addlistener(hl, 'ObjectBeingDestroyed', @(h, ev) lcl_obd_listener(hax2));

set(ancestor(hax2, 'figure'), 'CurrentAxes', hax2);

hl = line(w, angle(h), 'Parent', hax2);

ylabel(hax2, getString(message('signal:dspdata:dspdata:Phaseradians')));

set(ancestor(hax2, 'figure'), 'CurrentAxes', hax);

set(hl, 'Color', getcolorfromindex(hax2, 2));

setcoincidentgrid([hax hax2]);

set([hax hax2], ...
    'XGrid', 'On', ...
    'YGrid', 'On', ...
    'Box',   'On', ...
    'XLim',  [min(w) max(w)]);

if nargout
    varargout = {h};
end

% -------------------------------------------------------------------------
function lcl_obd_listener(hax2)

if ishandle(hax2)
    delete(hax2);
end

% [EOF]
