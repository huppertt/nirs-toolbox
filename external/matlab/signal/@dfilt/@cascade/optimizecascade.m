function Href = optimizecascade(this,Href,fn,varargin)
% This should be a private method.

%   Copyright 2009 The MathWorks, Inc.


% Test if Fixed-Point Designer is installed
if ~isfixptinstalled,
    error(message('signal:dfilt:cascade:optimizecascade:fixptTbxRq'));
end

% Test if response type is supported
iscoeffwloptimizable(this);

N = nstages(Href);
if iscell(fn),
    WL = fn{2};
    fn = fn{1};
    
    if length(WL) == 1,
        WL = WL*ones(N,1);
    elseif length(WL) ~= N,
        error(message('signal:dfilt:cascade:optimizecascade:wlVecWrongSize'));
    end
    try
        for i=1:N,
            Href.Stage(i) = fn(Href.Stage(i),WL(i),varargin{:});
        end
    catch ME,
        throwAsCaller(ME);
    end
else
    try
        for i=1:N,
            Href.Stage(i) = fn(Href.Stage(i),varargin{:});
        end
    catch ME,
        throwAsCaller(ME);
    end
end

setfdesign(Href,getfdesign(this));
setfmethod(Href, getfmethod(this));
