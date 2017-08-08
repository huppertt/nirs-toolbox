function this = notchbwapast(varargin)
%NOTCHBWAPAST   Construct a NOTCHBWAPAST object.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

this = fspecs.notchbwapast;

set(this, 'ResponseType', 'Notching Filter');

this.setspecs(varargin{:});


% [EOF]
