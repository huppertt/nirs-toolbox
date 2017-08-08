function [F, A] = getmask(this, fcns, ~, specs)
%GETMASK Get the mask.

%   Copyright 2005-2011 The MathWorks, Inc.

% If the specs were not passed in or are [], use the design specifications.
if nargin < 4 || isempty(specs)
    specs = getspecs(this.CurrentSpecs);
end

fpass = specs.Fpass;
fstop = specs.Fstop;

% Calculate the inverse sinc frequency vector and amplitude
isincF = linspace(0, fpass, 100);
isincA = 1./sinc(isincF*specs.FrequencyFactor).^specs.Power;

% Frequency vector is similar to Lowpass, but with the isincF inserted.
F = [1 fpass fliplr(isincF) isincF fstop fstop 1]*fcns.getfs()/2;

% Get the formatted Apass and Astop
apass = fcns.formatapass(specs.Apass);
astop = fcns.formatastop(specs.Astop);

% Convert the isinc amplitude to the correct units.
switch lower(fcns.getunits())
    case 'db'
        isincA = 20*log10(isincA);
    case 'squared'
        isincA = isincA.^2-1;
    case {'linear', 'zerophase'}
        isincA = isincA-1;
end

% Construct and amplitude vector.
A = [astop(1:2) apass(1)+fliplr(isincA) apass(2)+isincA ...
    apass(2)+isincA(end) astop(3:4)];

% [EOF]
