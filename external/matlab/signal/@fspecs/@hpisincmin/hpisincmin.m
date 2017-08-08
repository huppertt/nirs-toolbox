function this = hpisincmin(varargin)
%HPISINCMIN Construct a HPISINCMIN object.

%   Copyright 2011 The MathWorks, Inc.

% Override factory defaults inherited from lowpass
if nargin < 1,
    varargin{1} = .45;
end
if nargin < 2,
    varargin{2} = .55;
end

this = fspecs.hpisincmin;

respstr = 'Inverse-sinc highpass';
fstart = 1;
fstop = 2;
nargsnoFs = 4;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
