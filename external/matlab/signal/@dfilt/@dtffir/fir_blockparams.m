function s = fir_blockparams(Hd, mapstates,forceDigitalFilterBlock)
%FIR_BLOCKPARAMS   Return the fir specific block parameters.

%   Copyright 1988-2012 The MathWorks, Inc.

if nargin < 3
  forceDigitalFilterBlock = false;
end

s = blockparams(Hd.filterquantizer,forceDigitalFilterBlock);

if ~forceDigitalFilterBlock
    
  s.Coefficients = mat2str(get(reffilter(Hd), 'Numerator'),18);
  
  if strcmpi(mapstates, 'on')
    ic = getinitialconditions(Hd);
    if isempty(ic)
      ic = 0;
    end
    s.InitialStates = mat2str(ic);
  end
  
else
  
  % Target a Digital Filter block
  s.TypePopup = 'FIR (all zeros)';
  
  s.NumCoeffs = mat2str(get(reffilter(Hd), 'Numerator'),18);
  
  % IC
  if strcmpi(mapstates, 'on'),
    s.IC = mat2str(getinitialconditions(Hd));
  end  
end