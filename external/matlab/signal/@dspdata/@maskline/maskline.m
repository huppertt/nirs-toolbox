function this = maskline(varargin)
%MASKLINE   Construct a MASKLINE object.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

this = dspdata.maskline;

if nargin
    set(this, varargin{:});
end

% [EOF]
