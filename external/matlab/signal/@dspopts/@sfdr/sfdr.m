function this = sfdr(varargin)
%SFDR   Construct a SFDR options object

%   Copyright 2007 The MathWorks, Inc.

this = dspopts.sfdr;

if nargin   
    set(this, varargin{:});
end
% [EOF]
