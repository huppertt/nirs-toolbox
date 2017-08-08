function s = blockparams(Hd, mapstates, forceDigitalFilterBlock)
%BLOCKPARAMS Returns the parameters for BLOCK

% Copyright 1988-2012 The MathWorks, Inc.

if nargin < 3
  forceDigitalFilterBlock = false;
end

s = iir_blockparams(Hd,forceDigitalFilterBlock);

if ~forceDigitalFilterBlock
  
  s.FilterStructure = 'Direct Form II transposed';
  
  if strcmpi(mapstates, 'on')
    ic = getinitialconditions(Hd);
    if isempty(ic)
      ic = 0;
    end
    s.InitialStates = mat2str(ic);
  end
  
else
  
  s.IIRFiltStruct = 'Direct Form II transposed';
  
  if strcmpi(mapstates, 'on'),
    s.IC = mat2str(getinitialconditions(Hd));
  end
  
end