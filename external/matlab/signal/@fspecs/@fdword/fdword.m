function this = fdword(varargin)
%FDWORDER   Construct a FDWORDER object.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

this = fspecs.fdword;
this.ResponseType = 'Fractional Delay with Filter Order';
this.FilterOrder = 3;
if nargin>0,
    this.FilterOrder = varargin{1};
end

% [EOF]
