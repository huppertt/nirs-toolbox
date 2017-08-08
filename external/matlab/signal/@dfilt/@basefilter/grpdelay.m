function varargout = grpdelay(Hb,varargin)
%GRPDELAY Group delay of a discrete-time filter.

%   Copyright 1988-2004 The MathWorks, Inc.

if nargout,
    [Gd, w] = base_resp(Hb, 'computegrpdelay', varargin{:});
    varargout = {Gd, w};
else,
    [Hb, opts] = freqzparse(Hb, varargin{:});
    fvtool(Hb, 'grpdelay', opts);    
end

% [EOF]
