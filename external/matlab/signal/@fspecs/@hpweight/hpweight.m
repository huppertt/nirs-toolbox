function this = hpweight(varargin)
%HPWEIGHT   Construct a HPWEIGHT object.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

this = fspecs.hpweight;

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

respstr = 'Highpass';
fstart = 2;
fstop = 3;
nargsnoFs = 5;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
