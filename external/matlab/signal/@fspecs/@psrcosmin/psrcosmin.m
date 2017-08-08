function this = psrcosmin(varargin)
%PSRCOSMIN Construct a PSRCOSMIN object

%   Copyright 2008 The MathWorks, Inc.

this = fspecs.psrcosmin;

this.ResponseType = 'Minimum order raised cosine pulse shaping';

this.Astop = 60;

this.setspecs(varargin{:});

% [EOF]
