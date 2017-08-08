function realizemdl(Hd,varargin)
%REALIZEMDL Filter realization (Simulink diagram).
%     REALIZEMDL(Hd) automatically generates architecture model of filter
%     Hd in a Simulink subsystem block using individual sum, gain, and
%     delay blocks, according to user-defined specifications.
%
%     REALIZEMDL(Hd, PARAMETER1, VALUE1, PARAMETER2, VALUE2, ...) generates
%     the model with parameter/value pairs.
%
%     See also DFILT/REALIZEMDL

%   Copyright 2009-2011 The MathWorks, Inc.

grouper = SLM3I.ScopedGroupTransactions(); %#ok<NASGU>

% Parse input
[hTar,doMapCoeffsToPorts] = local_parseinput(Hd,varargin{:});

% Create model
pos = createmodel(hTar); 

% Generate filter architecture
DGDF = dgdfgen(Hd,hTar,doMapCoeffsToPorts);   % method defined in each DFILT class

% Expand dg_dfilt structure into directed graph
DG = expandToDG(DGDF,doMapCoeffsToPorts);

% Optimize direct graph
DG = optimizedg(Hd,hTar,DG);

% Generate mdl system
dg2mdl(DG,hTar,pos);

% Refresh connections
refreshconnections(hTar);

% Open system
opengeneratedmdl(hTar);

%------------------------------------------------------
function [hTar,doMapCoeffsToPorts] = local_parseinput(Hd,varargin)

% check if filter realizable
if ~isrealizable(Hd)
     error(message('signal:dfilt:abstractfilter:realizemdl:Notsupported', Hd.FilterStructure));
end

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
    error(message('signal:dfilt:abstractfilter:realizemdl:InvalidParameter'));
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
