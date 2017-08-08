function scaleflag = determinescaleflag(d)
%DETERMINESCALEFLAG Determine the correct scaling flag.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Determine whether to scale passband or not
switch get(d,'PassbandScale'),
case 'on',
    scaleflag = 'scale';
case 'off'
    scaleflag = 'noscale';
end
