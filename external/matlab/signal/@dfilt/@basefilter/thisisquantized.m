function b = thisisquantized(this)
%THISISQUANTIZED   Returns true if the filter is not set to double.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

% Get the property handle so that we can check to see if FilterDesign is
% installed.  If it isn't just assume we aren't quantized.
p = findprop(this, 'Arithmetic');
if isempty(p)
    b = false;
else
    if strcmpi(p.AccessFlags.PublicGet, 'Off')
        b = false;
    else
        b = ~strcmpi(this.Arithmetic, 'double');
    end
end

% [EOF]
