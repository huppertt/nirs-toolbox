function specs = whichspecs(h)
%WHICHSPECS

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

specs = ft_whichspecs(h);

specs(end+1).name   = 'invSincFreqFactor';
specs(end).datatype = 'udouble';
specs(end).defval   = 1;
specs(end).callback = [];
specs(end).descript = 'magspec';

specs(end+1).name   = 'invSincPower';
specs(end).datatype = 'udouble';
specs(end).defval   = 1;
specs(end).callback = [];
specs(end).descript = 'magspec';

% [EOF]
