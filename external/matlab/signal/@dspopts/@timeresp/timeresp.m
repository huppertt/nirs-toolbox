function this = timeresp(varargin)
%TIMERESP   Construct a TIMERESP object.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

this = dspopts.timeresp;

if nargin
    set(this, varargin{:});
end

% [EOF]
