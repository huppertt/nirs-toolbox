function [F, A] = getmask(this, fcns, rcf, specs)
%GETMASK   Get the mask.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

% If the specs were not passed in or are [], use the design specifications.
if nargin < 4 || isempty(specs)
    specs = getspecs(this.CurrentSpecs);
    fpass = specs.Fpass;
    fstop = specs.Fstop;
else
    fpass = specs.Fpass;
    fstop = specs.Fstop;
    if ~specs.NormalizedFrequency
        fpass = fpass/specs.Fs*2;
        fstop = fstop/specs.Fs*2;
    end
end

% The frequency vector is always the same.
F = [0 fpass fpass 1 1 fstop fstop 0]*fcns.getfs()/2;

A = fcns.gethighlow(specs);

% [EOF]
