function Hd = ellip(this, varargin)
%ELLIP Elliptic digital filter design.

%   Copyright 2007 The MathWorks, Inc.

Hd = design(this, 'ellip', varargin{:});
h = getfmethod(Hd);

if ishp(this),
    if isa(Hd,'mfilt.abstractmultirate'),
        error(message('signal:fspecs:hbmin:ellip:InvalidStructure'));
    end
    Hd = iirlp2hp(Hd,.5,.5);
    % Reset the contained FMETHOD.
    Hd.setfmethod(h);
end

% [EOF]
