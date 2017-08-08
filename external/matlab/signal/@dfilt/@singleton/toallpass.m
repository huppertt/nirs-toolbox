function h = toallpass(this)
%TOALLPASS   

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

if isallpass(this),
    [b,a] = tf(this);
    h = dfilt.allpass(a(2:end));
else
       error(message('signal:dfilt:singleton:toallpass:notAllpass'));
end

% [EOF]
