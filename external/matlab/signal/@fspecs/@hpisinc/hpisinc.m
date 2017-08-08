function this = hpisinc(varargin)
%HPISINC Construct a HPISINC object.

%   Copyright 2011 The MathWorks, Inc.

this = fspecs.hpisinc;

% Override factory defaults inherited from lowpass
if nargin < 3
    varargin{3} = .55;
    if nargin < 2
        varargin{2} = .45;
        if nargin < 1
            varargin{1} = 10;
        end
    end
end

respstr = 'Inverse-sinc Highpass';
fstart = 2;
fstop = 3;
nargsnoFs = 5;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
