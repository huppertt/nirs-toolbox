function setspecs(this, varargin)
%SETSPECS   Set the specs.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

% If Astop2 not specified, make equal to Astop1 if the later is specified
if length(varargin) > 4 && length(varargin) < 7,
    % If Apass not specified set default
    if length(varargin) < 6,
        varargin{6} = this.Apass;
    end
    varargin{7} = varargin{5};
end

aswfs_setspecs(this,varargin{:});

% [EOF]
