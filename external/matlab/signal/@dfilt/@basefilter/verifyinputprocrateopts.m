function verifyinputprocrateopts(Hd, hTar, type)
%VERIFYINPUTPROCRATEOPTS

%   Copyright 2012 The MathWorks, Inc.

allowedInputProcessing = {'columnsaschannels','elementsaschannels'};
allowedRateOptsForColumnasaschannels = {'enforcesinglerate','allowmultirate'};
allowedRateOptsForElementsaschannels = {'enforcesinglerate','allowmultirate'};

isMultirate = true;
if ~isa(Hd,'dfilt.cascade') && ~isa(Hd,'dfilt.parallel') 
  Hd1.Stage(1) = Hd;
  numStages = 1;
  isMultirate = ~all(Hd.privRateChangeFactor==[1 1]);       
else
  Hd1 = Hd;
  numStages = nstages(Hd1);
end

for i=1:numStages
  % allowedInputProcessing will never be empty, the other two
  % (allowedRateOptsForColumnasaschannels, and
  % allowedRateOptsForElementsaschannels) might be empty depending on the
  % filter. 
  [inputProc, rateOptionsFrameBased, rateOptionsSampleBased] = validinputprocrateoptions(Hd1.Stage(i),type);
  
  [~,idxA,~] = intersect(allowedInputProcessing,inputProc);
  allowedInputProcessing = allowedInputProcessing(sort(idxA));
  
  [~,idxA,~] = intersect(allowedRateOptsForColumnasaschannels,rateOptionsFrameBased);
  allowedRateOptsForColumnasaschannels = allowedRateOptsForColumnasaschannels(sort(idxA));
  
  [~,idxA,~] = intersect(allowedRateOptsForElementsaschannels,rateOptionsSampleBased);  
  allowedRateOptsForElementsaschannels = allowedRateOptsForElementsaschannels(sort(idxA));
end

if isempty(allowedRateOptsForColumnasaschannels) && isempty(allowedRateOptsForElementsaschannels)
  isMultirate = false;  
end

%------------------------------------------------------------------------
if ~hTar.IsInputProcessingSpecified && ~hTar.IsRateOptionSpecified  
  % User does not specify input processing or rate option
  if any(strcmpi(allowedInputProcessing,'columnsaschannels'))
    % If columns as channels is allowed, set input processing to colums as
    % channels. Set the first available rate option if any is available. 
    hTar.InputProcessing = 'columnsaschannels';
    if ~isempty(allowedRateOptsForColumnasaschannels)
      hTar.RateOption = allowedRateOptsForColumnasaschannels{1};
    end
  else
    % If elements as channels is allowed, set input processing to elements
    % as channels. Set the first available rate option if any is available.
    hTar.InputProcessing = 'elementsaschannels';
    if ~isempty(allowedRateOptsForElementsaschannels)
      hTar.RateOption = allowedRateOptsForElementsaschannels{1};  
    end
  end
%------------------------------------------------------------------------
elseif hTar.IsInputProcessingSpecified && ~hTar.IsRateOptionSpecified
  % User specifies input processing and does not specify rate option
  if ~strcmpi(hTar.InputProcessing,'inherited') && ~any(strcmpi(allowedInputProcessing,hTar.InputProcessing))
    % The specified input processing is not inherited and is not allowed so
    % we warn. Set input processing to the first valid value. 
    warning(message('signal:dfilt:basefilter:verifyinputprocrateopts:InputProcNotSupported',allowedInputProcessing{1}))
    hTar.InputProcessing = allowedInputProcessing{1};
  end  
  if strcmpi(hTar.InputProcessing,'inherited')
    % Specified input processing is 'inherited'. Set the first valid rate
    % options if any.
    if ~isempty(allowedRateOptsForColumnasaschannels)
      hTar.RateOption = allowedRateOptsForColumnasaschannels{1};
    elseif ~isempty(allowedRateOptsForElementsaschannels)
      hTar.RateOption = allowedRateOptsForElementsaschannels{1};
    end
  elseif strcmpi(hTar.InputProcessing,'columnsaschannels') && ~isempty(allowedRateOptsForColumnasaschannels)
    % User specified columns as channels (or we set it to correct an
    % invalid value). Set the first valid rate option if any.
    hTar.RateOption = allowedRateOptsForColumnasaschannels{1};
  elseif strcmpi(hTar.InputProcessing,'elementsaschannels') && ~isempty(allowedRateOptsForElementsaschannels)
    % User specified elements as channels (or we set it to correct an
    % invalid value). Set the first valid rate option if any.    
    hTar.RateOption = allowedRateOptsForElementsaschannels{1};
  end
%------------------------------------------------------------------------
elseif  ~hTar.IsInputProcessingSpecified && hTar.IsRateOptionSpecified  
  % User does not specify input processing and specifies rate option.
  if any(strcmpi(allowedInputProcessing,'columnsaschannels')) && any(strcmpi(allowedRateOptsForColumnasaschannels,hTar.RateOption))
    % If the combination of columns as channels an the specified rate
    % option is allowed, then set input processing to columns as channels. 
    hTar.InputProcessing = 'columnsaschannels';
  elseif any(strcmpi(allowedInputProcessing,'elementsaschannels')) && any(strcmpi(allowedRateOptsForElementsaschannels,hTar.RateOption))
    % If the combination of elements as channels an the specified rate
    % option is allowed, then set input processing to elements as channels.     
    hTar.InputProcessing = 'elementsaschannels';
  elseif  any(strcmpi(allowedInputProcessing,'columnsaschannels')) && ~any(strcmpi(allowedRateOptsForColumnasaschannels,hTar.RateOption))
    % If columns as channels is allowed, and the specified rate option is
    % not allowed, then warn and set rate option to an allowed value. Set
    % input processing to colums as channels. Do not warn if the filter is
    % not multirate. 
    if isMultirate && ~isempty(allowedRateOptsForColumnasaschannels)
      warning(message('signal:dfilt:basefilter:verifyinputprocrateopts:RateOptionNotSupported2',allowedRateOptsForColumnasaschannels{1},'columnsaschannels'))
    end
    hTar.InputProcessing = 'columnsaschannels';
    if ~isempty(allowedRateOptsForColumnasaschannels)
      hTar.RateOption = allowedRateOptsForColumnasaschannels{1};
    end
  elseif any(strcmpi(allowedInputProcessing,'elementsaschannels')) && ~any(strcmpi(allowedRateOptsForElementsaschannels,hTar.RateOption))
    % If elements as channels is allowed, and the specified rate option is
    % not allowed, then warn and set rate option to an allowed value. Set
    % input processing to elements as channels. Do not warn if the filter
    % is not multirate.
    if isMultirate && ~isempty(allowedRateOptsForElementsaschannels)
      warning(message('signal:dfilt:basefilter:verifyinputprocrateopts:RateOptionNotSupported2',allowedRateOptsForElementsaschannels{1},'elementsaschannels'))
    end
    hTar.InputProcessing = 'elementsaschannels';
    if ~isempty(allowedRateOptsForElementsaschannels)
      hTar.RateOption = allowedRateOptsForElementsaschannels{1};
    end
  end
%------------------------------------------------------------------------  
elseif hTar.IsInputProcessingSpecified && hTar.IsRateOptionSpecified    
  % User specifies input processing and rate option
  if ~strcmpi(hTar.InputProcessing,'inherited') && ~any(strcmpi(allowedInputProcessing,hTar.InputProcessing))
    % The specified input processing is not inherited and is not allowed so
    % we warn. Set input processing to the first valid value. 
    warning(message('signal:dfilt:basefilter:verifyinputprocrateopts:InputProcNotSupported',allowedInputProcessing{1}))
    hTar.InputProcessing = allowedInputProcessing{1};
  end
  if strcmpi(hTar.InputProcessing,'columnsaschannels') && ~any(strcmpi(allowedRateOptsForColumnasaschannels,hTar.RateOption))
    % If the combination of columns as channels and the specified rate
    % option is not allowed warn about the invalid rate option and set it
    % to the firs valid value. Do not warn if the filter is not multirate.
    if isMultirate && ~isempty(allowedRateOptsForColumnasaschannels)
      warning(message('signal:dfilt:basefilter:verifyinputprocrateopts:RateOptionNotSupported',allowedRateOptsForColumnasaschannels{1}))
    end
    if ~isempty(allowedRateOptsForColumnasaschannels)
      hTar.RateOption = allowedRateOptsForColumnasaschannels{1};
    end
  elseif strcmpi(hTar.InputProcessing,'elementsaschannels') && ~any(strcmpi(allowedRateOptsForElementsaschannels,hTar.RateOption))
    % If the combination of elements as channels and the specified rate
    % option is not allowed warn about the invalid rate option and set it
    % to the firs valid value. Do not warn if the filter is not multirate.    
    if isMultirate && ~isempty(allowedRateOptsForElementsaschannels)
      warning(message('signal:dfilt:basefilter:verifyinputprocrateopts:RateOptionNotSupported',allowedRateOptsForElementsaschannels{1}))
    end
    if ~isempty(allowedRateOptsForElementsaschannels)
      hTar.RateOption = allowedRateOptsForElementsaschannels{1};
    end
  end                          
end
