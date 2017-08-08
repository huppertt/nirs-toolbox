function [xlim, ylim] = thispassbandzoom(this, fcns, Hd, ~)
%THISPASSBANDZOOM Returns the limits of the passband zoom.

%   Copyright 2005-2011 The MathWorks, Inc.

% Get the mask information from the specifications.
[f, a] = getmask(this, fcns);

% Get the limits from the mask.
xlim_spec = [f(3) f(4)];
ylim_spec = [a(10) a(4)];

% If there is no dynamic range in the specifications, try the measurements.
if ~isempty(Hd)
    m = measure(Hd);
    
    if isempty(m.Apass)
        m = measure(Hd, 'Fpass1', m.F3dB1, 'Fpass2', m.F3dB2);
    end
    
    m = copy(m);
    
    m.normalizefreq(true); %false, fcns.getfs());
    
    [f, a] = getmask(this, fcns, Hd.nominalgain, m); %measure(Hd));
    xlim_meas = [f(3) f(4)];
    ylim_meas = [a(10) a(4)];
    
    % If the measurements failed we might have nans.
    if ~any(isnan(xlim_meas))   
        if xlim_meas(1) > xlim_spec(1), xlim(1) = xlim_meas(1);
        else                           xlim(1) = xlim_spec(1); end

        if xlim_meas(2) < xlim_spec(2), xlim(2) = xlim_meas(2);
        else                           xlim(2) = xlim_spec(2); end
    else
        xlim = xlim_spec;
    end
    
    if ylim_meas(2) > ylim_spec(2), ylim = ylim_meas;
    else                           ylim = ylim_spec; end
else
    xlim = xlim_spec;
    ylim = ylim_spec;
end

% Calculate the dynamic range of the xlimits.
dr_xlim = diff(xlim);

% Add padding to the xlimits based on the dynamic range.
xlim = xlim + dr_xlim*[-1 1]/10;

% [EOF]
