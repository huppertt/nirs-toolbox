function varargout = stepz(Hb, N, varargin)
%STEPZ  Discrete-time filter step response.

%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(1,3,nargin,'struct'));

if nargout,
    
    % Calculate the maximum impzlength if one is not provided.
    if nargin < 2, N = max(impzlength(Hb)); end

    [y,t]     = base_resp(Hb, 'computestepz', N, varargin{:});
    varargout = {y,t};
else
    
    if nargin > 1, varargin = {N, varargin{:}}; end
    
    % Launch FVTool
    [Hb, opts] = timezparse(Hb, varargin{:});
    fvtool(Hb, 'step', opts);
end

% [EOF]
