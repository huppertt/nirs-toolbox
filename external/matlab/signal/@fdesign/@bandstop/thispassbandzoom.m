function [xlim, ylim] = thispassbandzoom(this, fcns, Hd, ~)
%THISPASSBANDZOOM Returns the limits of the passband zoom.

%   Copyright 2005-2011 The MathWorks, Inc.

% Get the mask information from the specifications.
[~, a] = getmask(this, fcns);

% Get the limits from the mask.
ylim_specified = [a(3) a(4); a(10) a(9)];

if ~isempty(Hd)
    
    m = measure(Hd);
    
    opts = {};
    if isempty(m.Apass1)
        opts = {'Fpass1', m.F3dB1};
    end
    
    if isempty(m.Apass2)
        opts = {opts{:}, 'Fpass2', m.F3dB2}; %#ok<CCAT>
    end
    
    if ~isempty(opts)
        m = measure(Hd, opts{:});
    end
    
    m = copy(m);
    m.normalizefreq(true);
    
    [~, a] = getmask(this, fcns, Hd.nominalgain, m);
    
    ylim_measured = [a(3) a(4); a(10) a(9)];
    
    % Use whichever value is larger.
    if ylim_measured(1,2) > ylim_specified(1,2), ylim(1,:) = ylim_measured(1,:);
    else                                        ylim(1,:) = ylim_specified(1,:); end

    if ylim_measured(2,2) > ylim_specified(2,2), ylim(2,:) = ylim_measured(2,:); %#ok<*NASGU>
    else                                        ylim(2,:) = ylim_specified(2,:); end

    ylim = ylim_measured;
    
else
    ylim = ylim_specified;
end

xlim = [0 fcns.getfs()/2];

% We have 2 passbands so we take the bigger.
ylim = [min(ylim(:, 1)) max(ylim(:, 2))];

% [EOF]
