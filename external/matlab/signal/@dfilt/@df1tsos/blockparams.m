function s = blockparams(Hd, mapstates, varargin)
%BLOCKPARAMS Returns the parameters for BLOCK

% Copyright 1988-2012 The MathWorks, Inc.

s = super_blockparams(Hd);
s.IIRFiltStruct = 'Direct form I transposed';

% IC
if strcmpi(mapstates, 'on'),
   
    ic  = getinitialconditions(Hd);    
    s.ICNum = mat2str(ic.Num);
    s.ICDen = mat2str(ic.Den);
    
end
