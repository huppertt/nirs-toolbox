function [F, A] = getmask(this, fcns, rcf, specs)
%GETMASK   Get the mask.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

% If the specs were not passed in or are [], use the design specifications.
if nargin < 4 || isempty(specs)
    specs = getspecs(this.CurrentSpecs);
end

apass = specs.Apass;
fpass = specs.Fpass;
fstop = specs.Fstop;
astop = specs.Astop;

% Get low stopband points
stoplow = fcns.findbottom(-astop);

% The frequency vector is always the same.
F = [fpass fpass fstop fstop NaN fpass fstop]*fcns.getfs()/2;
 
% Get the ripples into the correct units.
apass = fcns.formatapass(apass);

if strcmpi(lower(fcns.getunits()), 'zerophase'), 
    apass = -apass; 
    stoplow = 0; 
end 

A = [stoplow apass(1) apass(1) stoplow NaN apass(2) apass(2)];

% [EOF]
