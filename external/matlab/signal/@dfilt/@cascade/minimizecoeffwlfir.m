function Hbest = minimizecoeffwlfir(this,Href,varargin)
%   This should be a private method.

%   Copyright 2009 The MathWorks, Inc.

Hbest = optimizecascade(this,Href,@minimizecoeffwl,varargin{:});

% [EOF]
