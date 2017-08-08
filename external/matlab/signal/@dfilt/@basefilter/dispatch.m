function Hd = dispatch(Hb)
%DISPATCH Returns the contained DFILT objects.

%   Author: V. Pellissier
%   Copyright 1988-2004 The MathWorks, Inc.

error(message('signal:dfilt:basefilter:dispatch:NotSupported', class( Hb )))

% [EOF]
