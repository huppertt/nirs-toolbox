function s = blockparams(Hd, mapstates, forceDigitalFilterBlock)
%BLOCKPARAMS Returns the parameters for BLOCK

% Copyright 1988-2012 The MathWorks, Inc.

if nargin < 3
  forceDigitalFilterBlock = false;
end

s = fir_blockparams(Hd, mapstates,forceDigitalFilterBlock);

if ~forceDigitalFilterBlock
  s.FilterStructure = 'Direct form symmetric';
else
  s.FIRFiltStruct = 'Direct form symmetric';
end

