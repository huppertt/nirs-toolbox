function s = blockparams(Hd, mapstates, varargin)
%BLOCKPARAMS   Return the block parameters.

%   Copyright 1988-2012 The MathWorks, Inc.

s = super_blockparams(Hd);
s.IIRFiltStruct = 'Direct form II';

% IC
if strcmpi(mapstates, 'on'),
    ic = getinitialconditions(Hd);
    s.IC  = mat2str(ic);
end
