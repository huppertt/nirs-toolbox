function Hd = kaiserwin(this, varargin)
%KAISERWIN   Design a kaiser-window filter.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

Hd = design(this, 'kaiserwin', varargin{:});


% [EOF]
