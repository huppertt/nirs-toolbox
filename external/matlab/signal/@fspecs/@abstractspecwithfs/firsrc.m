function varargout = firsrc(this, method, L, M, varargin)
%FIRSRC   Design an sample-rate converter.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

[struct varargin] = parsemultiratestruct(this, 'firsrc', varargin{:});

% Allow subclasses to overload this part of the process.  Some subclasses
% need to keep the length of the filter evenly divisible by the rate change
% factors.
[varargout{1:nargout}] = firmultirate(this, method, varargin{:});

% Get the coefficients from the Numerator of the returned object.
b = get(varargout{1}, 'Numerator');

% Create the multirate filter.
Hm = feval(['mfilt.' struct], L, M, L*b);

% Keep the FMETHOD for mask/measurement use later on.
setfmethod(Hm, getfmethod(varargout{1}));

varargout{1} = Hm;

% [EOF]
