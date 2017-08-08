function setspecs(this, varargin)
%SETSPECS   Set the specifications

%   Copyright 2008 The MathWorks, Inc.

% Find strings in varargin and set an error if a linear magnitude is
% specified. 
strs = varargin(cellfun(@ischar, varargin));

if ~isempty(intersect(strs,{'linear','squared'})),
    error(message('signal:fspecs:comblbwgbwnsh:setspecs:invalidSpecs'));
end

if nargin>1
    %combtype was prevously added in fdesign as the last argument of
    %varargin. Set CombType and delete it before setting the rest of the
    %specs
    this.CombType = varargin{end};
    varargin(end) = [];
end
aswfs_setspecs(this,varargin{:});

% [EOF]
