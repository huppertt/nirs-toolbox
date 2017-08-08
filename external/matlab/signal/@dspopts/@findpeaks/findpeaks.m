function this = findpeaks(varargin)
%FINDPEAKS Construct a FINDPEAKS options object

%   Copyright 2007 The MathWorks, Inc.

this = dspopts.findpeaks;

if nargin   
    set(this, varargin{:});
end

% [EOF]
