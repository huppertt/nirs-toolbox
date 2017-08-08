function constr = filtquant_plugins(h,arith)
%FILTQUANT_PLUGINS Table of filterquantizer plugins

%   Author(s): R. Losada
%   Copyright 1988-2005 The MathWorks, Inc.

switch arith
    case 'double',
        %#function dfilt.filterquantizer
        constr = 'dfilt.filterquantizer';
    case 'single',
        %#function dfilt.singlefilterquantizer
        constr = 'dfilt.singlefilterquantizer';
    otherwise,
        constr = thisfiltquant_plugins(h,arith);
end
