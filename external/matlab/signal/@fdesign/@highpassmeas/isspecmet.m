function b = isspecmet(this, hfdesign, args)
%ISSPECMET   True if the object is specmet.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

if nargin < 2
    hfdesign = get(this, 'Specification');
end

% Get the specifications from the FDesign object.
specs = measureinfo(hfdesign);

if nargin < 3,
    Apt  = 0;
    Ast  = 0;
else
    Apt  = args.Apasstol;
    Ast  = args.Astoptol;
    % Replace specs with those in args (these may have been slightly
    % modified due to the design method not actually generating a filter
    % that meets the specs)
    specs.Apass = args.Apass;
    specs.Astop = args.Astop;
end

% Return true if the measured Apass is less than or equal to the specificed
% and if the measured Astop is greater than or equal to the specified. Or
% if they are within tolerance specified.
b = (isempty(specs.Astop) || this.Astop >= specs.Astop || norm(specs.Astop-this.Astop,inf) <= Ast) ...
    && (isempty(specs.Apass) || this.Apass <= specs.Apass || norm(specs.Apass-this.Apass,inf) <= Apt);

% [EOF]
