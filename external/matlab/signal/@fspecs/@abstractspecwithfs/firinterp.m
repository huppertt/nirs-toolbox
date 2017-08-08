function varargout = firinterp(this, method, varargin)
%FIRINTERP   Create an FIRINTERP object.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

if nargin > 2 && isnumeric(varargin{1}) && isscalar(varargin{1})
    L = varargin{1};
    varargin(1) = [];
else
    L = getinterpfactor(this);
end

% Look for the structure string.  Default to 'firinterp'.
[struct varargin] = parsemultiratestruct(this, 'firinterp', varargin{:});

% FIRMULTIRATE takes care of 2L requirement for all multirate filters.
[varargout{1:nargout}] = firmultirate(this, method, varargin{:});

% Get the numerator and the interpolation factor which will be needed to
% built the multirate filter object.
b = get(varargout{1}, 'Numerator');

% Create the INTERP filter and make sure we multiply the coefficients by L.
Hm = feval(['mfilt.' struct], L, L*b);

setfmethod(Hm, getfmethod(varargout{1}));

varargout{1} = Hm;

% [EOF]
