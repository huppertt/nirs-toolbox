function [inputProc, rateOptionsFrameBased, rateOptionsSampleBased] = validinputprocrateoptions(H,type)
%VALIDINPUTPROCRATEOPTIONS
% Input TYPE can be 'block' or 'realizemdl'

%   Copyright 2012 The MathWorks, Inc.

% Get rules for input processing and rate options for the block and
% realizemdl methods. Structures that do not support one of the methods
% will return don't care values. The realizemdl and block methods should
% handle the checks for unsupported structures. 

% Start with all possible values for input processing and rate options.
inputProc = {'columnsaschannels','elementsaschannels'};
rateOptionsFrameBased = {'enforcesinglerate','allowmultirate'};
rateOptionsSampleBased = {'enforcesinglerate','allowmultirate'};

% Query the filters to see if what input processing and rate options they
% support.
if strcmpi(type,'realizemdl')
  
  % The doFrameProcessing method returns true if frame based processing is
  % supported in realizemdl for the filter at hand.
  if ~doFrameProcessing(H)
    inputProc(1) = [];
  end
  
  % The getrealizemdlraterestrictions method returns the rate options NOT
  % supported in realizemdl for the filter at hand.
  rFB = getrealizemdlraterestrictions(H,'columnsaschannels');
  rSB = getrealizemdlraterestrictions(H,'elementsaschannels');
  
else % block method
  
  % The getblockinputprocessingrestrictions method returns the input
  % processing options NOT supported in the block method for the filter at
  % hand.
  r = getblockinputprocessingrestrictions(H);
  if any(strcmp(r,'columnsaschannels'));
    idx = strcmp(inputProc,'columnsaschannels');  
    inputProc(idx) = [];
  end
  if any(strcmp(r,'elementsaschannels'));
    idx = strcmp(inputProc,'elementsaschannels');  
    inputProc(idx) = [];
  end
  
  % The getblockraterestrictions method returns the rate options NOT
  % supported in the block method for the filter at hand.
  rFB = getblockraterestrictions(H,'columnsaschannels');
  rSB = getblockraterestrictions(H,'elementsaschannels');
      
end

if any(strcmp(rFB,'enforcesinglerate'));
  idx = strcmp(rateOptionsFrameBased,'enforcesinglerate');
  rateOptionsFrameBased(idx) = [];
end
if any(strcmp(rFB,'allowmultirate'))
  idx = strcmp(rateOptionsFrameBased,'allowmultirate');
  rateOptionsFrameBased(idx) = [];
end

if any(strcmp(rSB,'enforcesinglerate'))
  idx = strcmp(rateOptionsSampleBased,'enforcesinglerate');
  rateOptionsSampleBased(idx) = [];
end
if any(strcmp(rSB,'allowmultirate'))
  idx = strcmp(rateOptionsSampleBased,'allowmultirate');
  rateOptionsSampleBased(idx) = [];
end

if ~any(strcmp(inputProc,'columnsaschannels'))
  rateOptionsFrameBased = {};
end

if ~any(strcmp(inputProc,'elementsaschannels'))
  rateOptionsSampleBased = {};
end
