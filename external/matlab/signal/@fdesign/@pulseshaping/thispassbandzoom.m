function [xlim, ylim] = thispassbandzoom(this, fcns, Hd, ~)
%PASSBANDZOOM   Returns the limits of the passband zoom.

%   Copyright 2011 The MathWorks, Inc.

% Get the mask information from the subclass.
[f, ~] = getmask(this, fcns);

% Compute y-axis zomm interval based on measured ripple
m =  measure(Hd);
Apass = m.Apass;
units = feval(fcns.getunits);

if strcmpi('linear',units)
  Apass = 10^(Apass/20);
  ylim = [2-Apass Apass];
elseif strcmpi('squared',units)
  Apass = 10^(Apass/10);
  ylim = [2-Apass Apass];
else % dB
  ylim = [-Apass Apass];
end

% Add a space of one tenth of the passband frequency to the x-axis zoom
% interval
xlim = [0 f(2)+f(2)/10];


% [EOF]
