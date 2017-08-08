function varargout = impz(Hb, N, varargin)
%IMPZ Impulse response of digital filter

%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(1,3,nargin,'struct'))

if nargout,
    if nargin < 2, N = max(impzlength(Hb)); end
   
    [y,t]     = base_resp(Hb, 'computeimpz', N, varargin{:});
    varargout = {y,t};
else
    if nargin > 1, varargin = {N, varargin{:}}; end

    [Hb, opts] = timezparse(Hb, varargin{:});    
    fvtool(Hb, 'impulse',opts);
end

% [EOF]
