function hopts = getOptshndl(h,arrayh)
%GETOPTSHNDL Get handle to frame with options.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Get handle to options frame
hopts = find(arrayh,'Tag','siggui.iirgrpdelayoptsframe');

