function [F, A] = getmask(this, fcns, ~, specs)
%GETMASK Get the mask.

%   Copyright 2011 The MathWorks, Inc.

% If the specs were not passed in or are [], use the design specifications.
if nargin < 4 || isempty(specs)
    specs = getspecs(this.CurrentSpecs);
end

fpass = specs.Fpass;
fstop = specs.Fstop;

% Calculate the inverse sinc frequency vector and amplitude
isincF = linspace(fpass, 1, 100);
isincA = (1./sinc((1-isincF)*specs.FrequencyFactor).^specs.Power);

% Frequency vector is similar to Highpass, but with the isincF inserted.
F = [0 fpass isincF fliplr(isincF) fstop fstop 0]*fcns.getfs()/2;

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
A = [astop(1:2)  apass(1)+isincA apass(2)+fliplr(isincA) ...
    apass(2)+isincA(1) astop(3:4)];

% [EOF]
