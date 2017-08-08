function Hd = optimizestopbandfir(this,Href,WL,varargin)
%OPTIMIZESTOPBANDFIR Optimize stopband.
%   This should be a private method.

%   Copyright 2009 The MathWorks, Inc.


Hd = optimizecascade(this,Href,{@maximizestopband,WL},varargin{:});


% [EOF]
