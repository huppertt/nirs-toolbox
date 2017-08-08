function Hd = butter(this, varargin)
%BUTTER   Butterworth digital filter design.

%   Copyright 1988-2005 The MathWorks, Inc.

Hd = design(this, 'butter', varargin{:});
h = getfmethod(Hd);

if ishp(this),
    if isa(Hd,'mfilt.abstractmultirate'),
        error(message('signal:fspecs:hbord:butter:InvalidStructure'));
    end
    Hd = iirlp2hp(Hd,.5,.5);
    % Reset the contained FMETHOD.
    Hd.setfmethod(h);
end
