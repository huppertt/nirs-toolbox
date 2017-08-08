function this = pssqrtrcosmin(varargin)
%PSRCOSMIN Construct a PSSQRTRCOSMIN object

%   Copyright 2008 The MathWorks, Inc.

this = fspecs.pssqrtrcosmin;

this.ResponseType = 'Minimum order square root raised cosine pulse shaping';

% This is the half of the default raised cosine stop band attenuation, which is
% 60 dB.
this.Astop = 30;

this.setspecs(varargin{:});

% [EOF]