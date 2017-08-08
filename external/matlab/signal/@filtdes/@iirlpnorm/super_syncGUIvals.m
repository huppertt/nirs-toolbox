function super_syncGUIvals(h, arrayh)
%SYNCGUIVALS Sync values from frames.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Call the super's method
lpnorm_syncGUIvals(h,arrayh);

% Call the num den order method
syncGUIvals(get(h,'numDenFilterOrderObj'),arrayh);

% Get handle to options frame
hopts = getOptshndl(h,arrayh);

set(h,'initNum',evaluatevars(get(hopts,'InitNum')));
set(h,'initDen',evaluatevars(get(hopts,'InitDen')));
