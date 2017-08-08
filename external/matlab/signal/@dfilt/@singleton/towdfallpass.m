function h = towdfallpass(this)
%TOWDFALLPASS   

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

if isallpass(this),
    [b,a] = tf(this);
    h = dfilt.wdfallpass(a(2:end));
else
       error(message('signal:dfilt:singleton:towdfallpass:notAllpass'));
end

% [EOF]
