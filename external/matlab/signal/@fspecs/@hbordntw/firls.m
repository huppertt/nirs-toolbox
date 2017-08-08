function Hd = firls(this, varargin)
%FIRLS   Design a least-squares FIR filter.
%   Hd = FIRLS(Hs) designs a least-squares FIR filter that meets the
%   specifications in Hs.

%   Author(s): R. Losada
%   Copyright 1988-2006 The MathWorks, Inc.

Hd = design(this, 'firls', varargin{:});
h = getfmethod(Hd);

if ishp(this),
    Hd = firlp2hp(Hd);
    % Reset the contained FMETHOD.
    Hd.setfmethod(h);
end

% [EOF]
