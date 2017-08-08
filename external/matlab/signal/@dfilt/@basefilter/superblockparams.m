function [lib, srcblk, s] = superblockparams(Hd, mapstates, link2obj, varname, hTar)
%SUPERBLOCKPARAMS

%   Copyright 2006-2012 The MathWorks, Inc.

try
  forceDigitalFilterBlock = false; 
  [lib, srcblk, hasInProc, hasRateOpts] = blocklib(Hd,link2obj);
  
  if hasInProc && strcmp(hTar.InputProcessing,'inherited') && ...
      any(strcmp({'Discrete FIR Filter','Discrete Filter','Allpole Filter'}, srcblk))
    % For inherited input processing force a DigitalFilter block. Set the
    % forceDigitalFilterBlock flag of the blocklib method to do this. 
    forceDigitalFilterBlock = true;
    [lib, srcblk, hasInProc, hasRateOpts] = blocklib(Hd,link2obj,forceDigitalFilterBlock);
  end  
  
catch %#ok<CTCH>
  error(message('signal:dfilt:basefilter:superblockparams:noBlock', varname));
end

if strcmpi(link2obj,'on'),
  s = objblockparams(Hd, varname);
else
  s = blockparams(Hd, mapstates,forceDigitalFilterBlock); 
end

if strcmp([lib '/' srcblk],'simulink/Discrete/Discrete Filter')
   s.a0EqualsOne = 'on';
end

verifyinputprocrateopts(Hd, hTar, 'block');

if hasInProc
  switch (hTar.InputProcessing)
    case 'columnsaschannels'
      InputProcessing = 'Columns as channels (frame based)';
    case 'elementsaschannels'
      InputProcessing = 'Elements as channels (sample based)';
    case 'inherited'
      InputProcessing = 'Inherited (this choice will be removed - see release notes)';
  end
  s.InputProcessing = InputProcessing;
end

if ~isspblksinstalled && any(strcmp({'inherited'},hTar.InputProcessing))
  error(message('signal:dfilt:basefilter:superblockparams:inputProcessingNotSupported'));
end

if hasRateOpts
  switch(hTar.RateOption)
    case 'enforcesinglerate'
      RateOption = 'Enforce single-rate processing';
    case 'allowmultirate'
      RateOption = 'Allow multirate processing';
  end
  s.(getRateOptionParamName(Hd)) = RateOption;
end
% [EOF]
