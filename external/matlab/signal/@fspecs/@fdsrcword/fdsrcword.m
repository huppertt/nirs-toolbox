function this = fdsrcword(varargin)
%FDSRCWORD Construct a FDSRCWORD object

%   Copyright 2007 The MathWorks, Inc.

this = fspecs.fdsrcword;
this.ResponseType = 'Farrow SRC with Polynomial Order';
this.PolynomialOrder = 3;
if nargin>0,
    this.PolynomialOrder = varargin{1};
end

% [EOF]
