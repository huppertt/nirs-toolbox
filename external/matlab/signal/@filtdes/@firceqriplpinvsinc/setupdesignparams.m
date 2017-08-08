function args = setupdesignparams(h, d)
%SETUPDESIGNPARAMS

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

args = firceqrip_setupdesignparams(h, d);

args = {args{:},'invsinc',...
        [get(d,'invSincFreqFactor'),get(d,'invSincPower')]};

% [EOF]
