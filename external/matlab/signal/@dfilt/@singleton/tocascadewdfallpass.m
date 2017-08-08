function h = tocascadewdfallpass(this)
%TOCASCADEWDFALLPASS   

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

if isallpass(this),
    [b,a] = tf(this);
    h = dfilt.cascadewdfallpass(a(2:end));
else
       error(message('signal:dfilt:singleton:tocascadewdfallpass:notAllpass'));
end

% [EOF]
