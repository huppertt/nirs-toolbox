function this = sosview(varargin)
%SOSVIEW   Construct a SOSVIEW object.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

this = dspopts.sosview;

if nargin
    set(this, varargin{:});
end

% [EOF]
