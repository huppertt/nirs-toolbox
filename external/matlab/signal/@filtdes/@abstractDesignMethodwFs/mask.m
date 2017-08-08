function varargout = mask(d, hax)
%MASK Draws the mask to an axes.

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

if nargin < 2, hax = gca; end

h = info2mask(maskinfo(d), hax);

if nargout, varargout = {h}; end

% [EOF]
