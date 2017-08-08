function [xlim, ylim] = thispassbandzoom(this, fcns, Hd, hfm) %#ok<INUSD>
%THISPASSBANDZOOM   Returns the limits of the passband zoom.

%   Copyright 2005-2014 The MathWorks, Inc.

% Get the mask information from the subclass.
[f, a] = getmask(this, fcns);

% Get the limits from the mask.
xlim_specified = [f(3) fcns.getfs()/2];
ylim_specified = [a(5) a(4)];

% Calculate the dynamic range of the ylimits.

% If there is no dynamic range in the specifications, try the measurements.
if ~isempty(Hd)
    
    m = measure(Hd);
    if isempty(m.Apass)
        m = measure(Hd, 'Fpass', m.F3dB);
    end
    
    [f, a] = getmask(this, fcns, Hd.nominalgain, m);
    xlim_measured = [f(3) fcns.getfs()/2];
    ylim_measured = [a(5) a(4)];

    if xlim_measured(1) > xlim_specified(1)
        xlim = xlim_measured;
    else
        xlim = xlim_specified;
    end
    
    % If the measured Apass is greater than that 
    if ylim_measured(2) > ylim_specified(2) || ...
            diff(ylim_specified) == 0 || any(isnan(ylim_specified))
        ylim = ylim_measured;
    else
        ylim = ylim_specified;
    end
    
else
    xlim = xlim_specified;
    ylim = ylim_specified;
end

% Calculate the dynamic range of the xlimits.
dr_xlim = diff(xlim);

% Add padding to the xlimits based on the dynamic range.
xlim(1) = xlim(1)-dr_xlim/10;

% [EOF]
