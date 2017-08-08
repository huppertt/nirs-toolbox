function setspecs(this, varargin)
%SETSPECS   Set the specs.

%   Copyright 2007 The MathWorks, Inc.

for i=2:nargin,
    % Set Type
    if strcmpi(varargin{i-1},'type'),
        this.Type = varargin{i};
        varargin(i-1:i) =[];
        break;
    end
end
abstract_setspecs(this, varargin{:});

% [EOF]
