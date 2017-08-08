function this = hpstopfpass(varargin)
%HPSTOPFPASS   Construct a HPSTOPFPASS object.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

% Override factory defaults inherited from lowpass
if nargin < 1,
    varargin{1} = 10;
end
if nargin < 2,
    varargin{2} = .45;
end
if nargin < 3,
    varargin{3} = .55;
end

this = fspecs.hpstopfpass;

respstr = 'Highpass with passband frequency';
fstart = 2;
fstop = 3;
nargsnoFs = 4;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
