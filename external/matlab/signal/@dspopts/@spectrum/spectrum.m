function this = spectrum(varargin)
%SPECTRUM   Options object for PSD and mean-square spectrum analysis.
%   
%   To create a PSD or mean-square spectrum options object use the spectrum
%   object methods <a href="matlab:help spectrum/psdopts">psdopts</a> and <a href="matlab:help spectrum/msspectrumopts">msspectrumopts</a>, respectively.
%
%   See also SPECTRUM, SPECTRUM/PSEUDOSPECTRUMOPTS

%   Author(s): J. Schickler 
%   Copyright 1988-2006 The MathWorks, Inc. 


this = dspopts.spectrum;

if nargin
    set(this, varargin{:});
end

% [EOF]
