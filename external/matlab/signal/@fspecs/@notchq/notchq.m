function this = notchq(varargin)
%NOTCHQ   Construct a NOTCHQ object.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

this = fspecs.notchq;

set(this, 'ResponseType', 'Notching Filter');

this.setspecs(varargin{:});

% [EOF]
