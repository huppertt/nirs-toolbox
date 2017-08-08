function varargout = firmultirate(this, method, varargin)
%FIRMULTIRATE   Design a multirate conforming to the 2L polyphase rule.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

[varargout{1:nargout}] = feval(method, this, varargin{:});

b = get(varargout{1}, 'Numerator');

L = get(this, 'Band');

% If the length of b is not divisible by 2L, create a specify order object
% and force the filter order to the next factor of 2L.
if rem(length(b),2*L)
    h  = fspecs.nyqordastop;
    set(h, ...
        'Band', L, ...
        'FilterOrder', 2*L*ceil(length(b)/(2*L)), ...
        'Astop', get(this, 'Astop'));

    [varargout{1:nargout}] = feval(method, h, varargin{:});
end

% [EOF]
