function this = notchqap(varargin)
%NOTCHQAP   Construct a NOTCHQAP object.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

this = fspecs.notchqap;

set(this, 'ResponseType', 'Notching Filter');

this.setspecs(varargin{:});


% [EOF]
