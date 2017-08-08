function this = notchbw(varargin)
%NOTCHBW   Construct a NOTCHBW object.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

this = fspecs.notchbw;

set(this, 'ResponseType', 'Notching Filter');

this.setspecs(varargin{:});

% [EOF]
