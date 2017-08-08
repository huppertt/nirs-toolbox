function checkrequiredlicense(Hd,hTar)
%CHECKREQUIREDLICENSE check required license for realizemdl

%   Copyright 2009 The MathWorks, Inc.

% Check if Simulink is installed
[b, ~, ~, errObj] = issimulinkinstalled;
if ~b
    error(errObj);
end

% [EOF]
