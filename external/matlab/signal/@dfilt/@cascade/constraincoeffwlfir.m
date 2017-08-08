function Hd = constraincoeffwlfir(this,Href,WL,varargin)
%CONSTRAINCOEFFWLFIR Constrain coefficient wordlength.
%   This should be a private method.

%   Copyright 2009 The MathWorks, Inc.


Hd = optimizecascade(this,Href,{@constraincoeffwl,WL},varargin{:});


% [EOF]
