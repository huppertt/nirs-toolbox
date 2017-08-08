function p = scalarblockparams(this)
%SCALARBLOCKPARAMS   Return the parameters for the gain block.

%   Author(s): J. Schickler
%   Copyright 1999-2011 The MathWorks, Inc.
%     

p.ParamDataTypeStr          = 'float(''single'')';
p.OutDataTypeStr            = 'float(''single'')';
p.LockScale                 = 'on' ;
p.RndMeth                   = 'Nearest';  % Zero|Nearest|Ceiling|Floor
p.DoSatur                   = 'on';

% [EOF]
