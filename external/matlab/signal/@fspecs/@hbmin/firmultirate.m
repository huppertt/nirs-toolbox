function varargout = firmultirate(this, method, varargin)
%FIRMULTIRATE   Design a multirate conforming to the 2L polyphase rule.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

% Design the minimum order filter.
[varargout{1:nargout}] = feval(method, this, varargin{:});

b = get(varargout{1}, 'Numerator');

% If the length of b is not divisible by 4, then create a specify order
% object and increment the order to the next order that is divisible by 4.
% This is done to ensure that the zeros fall on the correct spots.
if rem(length(b),4)
    h  = fspecs.hbordastop;
    set(h, ...
        'Astop', get(this, 'Astop'), ...
        'FilterOrder', 4*ceil(length(b)/4));

    [varargout{1:nargout}] = feval(method, h, varargin{:});
end

% [EOF]
