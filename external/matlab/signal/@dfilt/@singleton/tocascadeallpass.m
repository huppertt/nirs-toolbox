function h = tocascadeallpass(this)
%TOCASCADEALLPASS   

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

if isallpass(this),
    [b,a] = tf(this);
    h = dfilt.cascadeallpass(a(2:end));
else
       error(message('signal:dfilt:singleton:tocascadeallpass:notAllpass'));
end

% [EOF]
