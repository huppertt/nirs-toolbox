function [x, y] = passbandzoom(this, varargin)
%PASSBANDZOOM   Return the limits for the passbandzoom.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

x = zeros(length(this), 2);
y = x;

% Loop over all the filters and get each ones passbandzoom limits.
for indx = 1:length(this)
    
    hfdesign = getfdesign(this(indx));
    hfmethod = getfmethod(this(indx));
    
    if isempty(hfdesign) || isempty(hfmethod)
        x = [NaN NaN];
        y = [NaN NaN];
        return;
    else
        [x(indx, :), y(indx, :)] = passbandzoom(hfdesign, hfmethod, this(indx), varargin{:});
    end
end

% Calculate the minimum area we will need to show ALL the passband zooms.
x = [min(x(:, 1)) max(x(:, 2))];
y = [min(y(:, 1)) max(y(:, 2))];

% [EOF]
