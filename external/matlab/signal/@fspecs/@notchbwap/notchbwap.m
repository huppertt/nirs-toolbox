function this = notchbwap(varargin)
%NOTCHBWAP   Construct a NOTCHBWAP object.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

this = fspecs.notchbwap;

set(this, 'ResponseType', 'Notching Filter');

this.setspecs(varargin{:});


% [EOF]
