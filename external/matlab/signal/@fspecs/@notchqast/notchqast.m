function this = notchqast(varargin)
%NOTCHQAST   Construct a NOTCHQAST object.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

this = fspecs.notchqast;

set(this, 'ResponseType', 'Notching Filter');

this.setspecs(varargin{:});


% [EOF]
