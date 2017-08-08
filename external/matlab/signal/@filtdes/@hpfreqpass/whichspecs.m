function specs = whichspecs(h)
%WHICHSPECS Determine which specs are required for this class.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Call super's method
specs = fp1_whichspecs(h);

% Change default value to something above the Fs/4 so halfband highpass can
% use it
specs.defval = 14400;


