function this = pseudospectrum(varargin)
%PSEUDOSPECTRUM   Options object for pseudospectrum analysis.
%
%   To create a pseudospectrum options object use the spectrum object
%   method <a href="matlab:help spectrum/pseudospectrumopts">pseudospectrumopts</a>.
%
%   See also SPECTRUM, SPECTRUM/MSSPECTRUMOPTS, SPECTRUM/PSDOPTS.

%   Author(s): J. Schickler
%   Copyright 1988-2006 The MathWorks, Inc.

this = dspopts.pseudospectrum;

if nargin
    set(this, varargin{:});
end

% [EOF]
