function super_realizemdl_composite(Hd,varargin)
%SUPER_REALIZEMDL_COMPOSITE realize composite model

%   Copyright 2009-2011 The MathWorks, Inc.

% Parse input
[hTar,doMapCoeffsToPorts] = local_parseinput(Hd,varargin{:});

% Create model
pos = createmodel(hTar);

% Generate filter architecture
msg = dgdfgen(Hd,hTar,doMapCoeffsToPorts,pos);   %#ok<NASGU> % method defined in each DFILT class

% Refresh connections
refreshconnections(hTar);

% Optimisations
optimize_mdl(hTar, Hd);

% Open system
opengeneratedmdl(hTar);
    
% -------------------------------------------------------------
function optimize_mdl(hTar, Hd)

% Optimize zero gains
if strcmpi(hTar.OptimizeZeros, 'on'),
     optimizezerogains(hTar, Hd);
end

% Optimize unity gains
if strcmpi(hTar.OptimizeOnes, 'on'),
     optimizeonegains(hTar, Hd);
end

% Optimize -1 gains
if strcmpi(hTar.OptimizeNegOnes, 'on'),
     optimizenegonegains(hTar, Hd);
end

% Optimise delay chains
if strcmpi(hTar.OptimizeDelayChains, 'on'),
    optimizedelaychains(hTar);
end
%------------------------------------------------------
function [hTar,doMapCoeffsToPorts] = local_parseinput(Hd,varargin)

% Parse inputs to target
hTar = uddpvparse('dspfwiztargets.realizemdltarget', varargin{:});

% Check if input processing and rateoptions were specified. If the
% IsInputProcessingSpecified, and IsRateOptionSpecified properties of hTar
% are already set to true, it means that the InputProcessing and/or
% RateOption values have already been specified directly in the object. Do
% not change the status of these properties in this case.
if ~hTar.IsInputProcessingSpecified
  idx = find(strcmpi(varargin,'InputProcessing'), 1);
  if ~isempty(idx)
    hTar.IsInputProcessingSpecified = true;
  end
end
if ~hTar.IsRateOptionSpecified
    idx = find(strcmpi(varargin,'RateOption'), 1);
    if ~isempty(idx)
      hTar.IsRateOptionSpecified = true;
    end
end

% Verify valid input processing and rate options
verifyinputprocrateopts(Hd,hTar, 'realizemdl');

% Clear gains and delays
hTar.gains = [];
hTar.delays = [];

% Check if the required license is installed
checkrequiredlicense(Hd,hTar)

% Check MapCoeffsToPorts state and coefficient names
if ~strcmpi(hTar.MapCoeffsToPorts,'on')&&~isempty(hTar.CoeffNames),
    error(message('signal:dfilt:basefilter:super_realizemdl_composite:InvalidParameter'));
end

% Parse parameters to hTar
try
    % Get coefficient names and values
    [hTar,doMapCoeffsToPorts] = parse_coeffstoexport(Hd,hTar);
    
    % Set filter states to hTar
    hTar = parse_filterstates(Hd,hTar);
catch ME
    throwAsCaller(ME);
end


% [EOF]
