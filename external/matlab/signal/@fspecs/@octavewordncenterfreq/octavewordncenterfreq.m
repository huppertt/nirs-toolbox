function this = octavewordncenterfreq(varargin)
%OCTAVEWORDNCENTERFREQ   Construct an OCTAVEWORDNCENTERFREQ object.

%   Author(s): V. Pellissier
%   Copyright 2006 The MathWorks, Inc.

this = fspecs.octavewordncenterfreq;
this.FilterOrder = 6;
if nargin>0,
    this.FilterOrder = varargin{1};
end
if nargin>1,
    this.F0 = varargin{2};
end
% [EOF]
