function [nfft, unitcircle] = freqzparse(hObj, varargin)
%FREQZPARSE Returns the inputs for freqz

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

nfft       = 512;
unitcircle = 1;

if nargin > 2,
    unitcircle = strmatch(lower(varargin{2}), {'half', 'whole', 'fftshift'});
end
if nargin > 1,
    nfft = varargin{1};
end

if length(nfft) > 1,
    unitcircle = 4;
end

% [EOF]
