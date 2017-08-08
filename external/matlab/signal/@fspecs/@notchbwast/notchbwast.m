function this = notchbwast(varargin)
%NOTCHBWAST   Construct a NOTCHBWAST object.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

this = fspecs.notchbwast;

set(this, 'ResponseType', 'Notching Filter');

this.setspecs(varargin{:});

% [EOF]
