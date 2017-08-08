function lpnorm_syncGUIvals(h, arrayh)
%SYNCGUIVALS Sync values from frames.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Call the super's method
base_syncGUIvals(h,arrayh);

% Get handle to options frame
hopts = getOptshndl(h,arrayh);

set(h,'Pnorm',evaluatevars(get(hopts,'PNormEnd')));
set(h,'initPnorm',evaluatevars(get(hopts,'PNormStart')));
set(h,'DensityFactor',evaluatevars(get(hopts,'DensityFactor')));

