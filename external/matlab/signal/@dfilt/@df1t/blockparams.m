function s = blockparams(Hd, mapstates, forceDigitalFilterBlock)
%BLOCKPARAMS Returns the parameters for BLOCK

% Copyright 1988-2012 The MathWorks, Inc.

if nargin < 3
  forceDigitalFilterBlock = false;
end

s = iir_blockparams(Hd,forceDigitalFilterBlock);

if ~forceDigitalFilterBlock
  s.FilterStructure = 'Direct Form I transposed';
  
  if strcmpi(mapstates, 'on')
    ic = getinitialconditions(Hd);
    
    if isempty(ic.Num)
      ic.Num = 0;
    end
    s.InitialStates = mat2str(ic.Num);
    
    if isempty(ic.Den)
      ic.Den = 0;
    end
    s.InitialDenominatorStates = mat2str(ic.Den);
  end
  
else
  
  s.IIRFiltStruct = 'Direct Form I transposed';
  
  if strcmpi(mapstates, 'on'),
    ic = getinitialconditions(Hd);
    
    s.ICDen = mat2str(ic.Den);
    s.ICNum = mat2str(ic.Num);
  end
  
end