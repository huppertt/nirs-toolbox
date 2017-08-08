function this = notchqapast(varargin)
%NOTCHQAPAST   Construct a NOTCHQAPAST object.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

this = fspecs.notchqapast;

set(this, 'ResponseType', 'Notching Filter');

this.setspecs(varargin{:});

% [EOF]
