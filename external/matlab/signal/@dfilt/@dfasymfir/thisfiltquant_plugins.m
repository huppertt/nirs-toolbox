function constr = thisfiltquant_plugins(h,arith)
%FILTQUANT_PLUGINS Table of filterquantizer plugins

%   Author(s): R. Losada
%   Copyright 1988-2005 The MathWorks, Inc.

switch arith
    case 'fixed',
        [wstr wid] = lastwarn;
        w = warning('off');
        constr = 'quantum.fixeddfsymfirfilterq';
        lastwarn(wstr, wid)
        warning(w);
end
