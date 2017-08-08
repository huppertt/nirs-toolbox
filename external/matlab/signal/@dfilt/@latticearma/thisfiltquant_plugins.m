function constr = thisfiltquant_plugins(h,arith)
%FILTQUANT_PLUGINS Table of filterquantizer plugins

%   Author(s): V. Pellissier
%   Copyright 1999-2005 The MathWorks, Inc.

switch arith
    case 'fixed',
        constr = 'quantum.fixedlatticearmafilterq';
end
