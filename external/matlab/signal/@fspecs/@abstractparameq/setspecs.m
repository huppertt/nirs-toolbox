function setspecs(this, varargin)
%SETSPECS   Set the specifications

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

% Find strings in varargin
strs = varargin(cellfun(@ischar, varargin));

if ~isempty(intersect(strs,{'linear','squared'})),
    error(message('signal:fspecs:abstractparameq:setspecs:invalidSpecs'));
end
aswfs_setspecs(this,varargin{:});

% [EOF]
