function s = blockparams(Hd, mapstates, forceDigitalFilterBlock)
%BLOCKPARAMS Returns the parameters for BLOCK

% Copyright 1988-2012 The MathWorks, Inc.

if nargin < 3
  forceDigitalFilterBlock = false;
end

s = blockparams(Hd.filterquantizer,forceDigitalFilterBlock);

if ~forceDigitalFilterBlock
  s.FilterStructure = 'Lattice AR';
  
  refHd = reffilter(Hd);
  
  s.Coefficients = mat2str(refHd.Lattice,18);
  
  if strcmpi(mapstates, 'on')
    ic = getinitialconditions(Hd);
    if isempty(ic)
      ic = 0;
    end
    s.InitialStates = mat2str(ic);
  end
  
else
  
  s.TypePopup = 'IIR (all poles)';
  s.AllPoleFiltStruct = 'Lattice AR';
  
  s.LatticeCoeffs = mat2str(Hd.Lattice,18);
  
  if strcmpi(mapstates, 'on'),
    s.IC = mat2str(getinitialconditions(Hd));
  end
  
end