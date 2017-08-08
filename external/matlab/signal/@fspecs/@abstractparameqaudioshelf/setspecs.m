function setspecs(this, varargin)
%SETSPECS   Set the specifications

%   Copyright 2008 The MathWorks, Inc.

%Set default values (different from the default fdesign.parameq equalizer
%values). If we don't do this, other default values are written into the
%object later on.
if nargin < 2
    this.FilterOrder = 2;
    this.F0 = 0;
    this.G0 = 10;
end

% Find strings in varargin
strs = varargin(cellfun(@ischar, varargin));

if ~isempty(intersect(strs,{'linear','squared'})),
    error(message('signal:fspecs:abstractparameqaudioshelf:setspecs:invalidSpecs'));
end

aswfs_setspecs(this,varargin{:});

% [EOF]
