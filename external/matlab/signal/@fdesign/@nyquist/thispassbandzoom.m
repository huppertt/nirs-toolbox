function [xlim, ylim] = thispassbandzoom(this, fcns, Hd, hfm)
%THISPASSBANDZOOM   

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

% Get the mask information from the subclass.
[f, a] = getmask(this, fcns);

% Get the limits from the mask.
xlim_specified = [0 f(3)];
ylim_specified = [a(5) a(4)];

% If there is no dynamic range in the specifications, try the measurements.
if ~isempty(Hd)
    m = measure(Hd);
    if isempty(m.Fpass)
        m = measure(Hd, 'Fpass', m.F3dB);
    end
    [f, a] = getmask(this, fcns, Hd.nominalgain, m);
    xlim_measured = [0 f(3)];
    ylim_measured = [a(5) a(4)];

    % Use the first frequency point.  Sometimes xlim_specified does not
    % contain Fpass because it was not actually specified.
    if xlim_measured(2) < xlim_specified(2)
        xlim = xlim_measured;
    else
        xlim = xlim_specified;
    end
    
    % If the specified ylim gives us the same two numbers, no ylim was
    % actually specified, use the Measured.  If the measured ylim is less
    % than specified we want to use the measured too.
    if diff(ylim_measured) ~= 0
        ylim = ylim_measured;
    elseif diff(ylim_specified) ~= 0
        ylim = ylim_specified;
    else
        % If neither of the ylims are usable return nans so we know to
        % disable.
        ylim = [NaN NaN];
    end
else
    xlim = xlim_specified;
    ylim = ylim_specified;
end

% Calculate the dynamic range of the limits.
dr_xlim = diff(xlim);

% Add padding to the xlimits based on the dynamic range.
xlim(2) = xlim(2)+dr_xlim/10;

% [EOF]
