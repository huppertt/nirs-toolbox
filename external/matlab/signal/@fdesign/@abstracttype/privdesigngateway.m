function varargout = privdesigngateway(this, method, varargin)
%PRIVDESIGNGATEWAY   Gateway for all of the design methods.
%   PRIVDESIGNGATEWAY(H, METHOD, PV-pairs)

%   Copyright 2005-2011 The MathWorks, Inc.

% We assume here that "method" is valid.

% Force a single output from THISDESIGN.
n = nargout;

if isdesignmethod(this, method)
    [varargout{1:n}] = thisdesign(this, method, varargin{:});
else
    error(message('signal:fdesign:abstracttype:privdesigngateway:invalidMethod', upper( method ), this.Specification));
end

% Store fdesign obj in dfilt
setfdesign(varargout{1},this);

% Store the design method string
setdesignmethod(varargout{1},method);

% Convert to System object if it has been requested
fm = getfmethod(varargout{1});
if nargout && isprop(fm,'SystemObject') && fm.SystemObject
  varargout{1} = sysobj(varargout{1});
end

if ~nargout
    Hd = varargout{1};
    varargout = {};
    if this.NormalizedFrequency,
        inputs = {'NormalizedFrequency', 'On'};
    else
        inputs = {'Fs', this.Fs};
    end
    inputs = {inputs{:}, 'DesignMask', 'on'};
    fvtool(Hd, inputs{:});
end

% [EOF]
