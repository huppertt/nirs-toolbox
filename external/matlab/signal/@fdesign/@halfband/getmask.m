function [F, A] = getmask(this, fcns, rcf, specs)
%GETMASK   Get the mask.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

% If the specs were not passed in or are [], use the design specifications.
if nargin < 4 || isempty(specs)
    specs = getspecs(this.CurrentSpecs);
end

fpass = specs.Fpass;
fstop = specs.Fstop;

% The frequency vector is always the same.
if strcmpi(this.Type,'Lowpass'),
    F = [1 fpass fpass 0 0 fstop fstop 1]*fcns.getfs()/2;
elseif strcmpi(this.Type,'Highpass'),
    F = [0 fpass fpass 1 1 fstop fstop 0]*fcns.getfs()/2;
end

A = fcns.gethighlow(specs);

% [EOF]
