function hm = measure(this, Hd, varargin)
%MEASURE   Measure this object.

%   Author(s): J. Schickler
%   Copyright 2005-2010 The MathWorks, Inc.

% Measure is only available to DSP System Toolbox license holders
supercheckoutfdtbxlicense(this)

error(nargchk(2,inf,nargin,'struct'));

if isempty(getfdesign(Hd))
    Hd = copy(Hd);
    setfdesign(Hd,this);
end

hm = feval(getmeasureconstructor(this), Hd, this, varargin{:});

% [EOF]
