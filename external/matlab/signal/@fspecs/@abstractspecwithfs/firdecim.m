function varargout = firdecim(this, method, varargin)
%FIRDECIM   Design a Decimation filter.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

% The window method might send in a window vector as the first input.  This
% is the only exception.
if nargin > 2 && isnumeric(varargin{1}) && isscalar(varargin{1})
    M = varargin{1};
    varargin(1) = [];
else
    M = getdecimfactor(this);
end

% Look for the structure string.  Default to 'firdecim'.
[struct varargin] = parsemultiratestruct(this, 'firdecim', varargin{:});

% FIRMULTIRATE takes care of 2L requirement for all multirate filters.
[varargout{1:nargout}] = firmultirate(this, method, varargin{:});

% Get the numerator and the interpolation factor which will be needed to
% built the multirate filter object.
b = get(varargout{1}, 'Numerator');

% Create the DECIM filter.
Hm = feval(['mfilt.' struct], M, b);

setfmethod(Hm, getfmethod(varargout{1}));

varargout{1} = Hm;

% [EOF]
