function Hs = tosysobj(this,returnSysObj)
%TOSYSOBJ Convert to a System object

%   Copyright 2013-2014 The MathWorks, Inc.

if ~returnSysObj
  % If returnSysObj is false, then it means that we want to know if the
  % System object conversion is supported for the class at hand. 
  
  Hs = true;
  return;
end  

t =  getfmethod(this);
coupledAllPassDesign = false;
if(isprop(t, 'FilterStructure'))
    desString = t.FilterStructure;
    switch(desString)
        case 'cascadeallpass',
            Hs = dsp.CoupledAllpassFilter('Minimum multiplier');
            coupledAllPassDesign = true;
        case 'cascadewdfallpass',
            Hs = dsp.CoupledAllpassFilter('Wave Digital Filter');
            coupledAllPassDesign = true;
    end
end

if coupledAllPassDesign
    % Sets properties of dsp.CoupledAllpassFilter matching filter structure in
    % source dfilt.cascade, including
    % - Coefficients (for branches with dfilt.cascade(wdf)allpass
    % - Presence of delays in branches, with various implementations in source
    %   dfilt.cascade, including through dfilt.delay and dfilt.dffir
    % - Individual branch gains, with various implementations in source
    %   dfilt.cascade, including through dfilt.scalar and dfilt.dffir
    setProperties(this, Hs);
    
    % Generates a warning if source dfilt has persistent memory and non-trivial
    % internal states, since dsp.CoupledAllpassFilter does not support nonzero
    % initial conditions
    warnIfStatesNontrivial(this);
else
    % Create a cascade of System objects
    if this.nstages == 0
        error(message('signal:dfilt:cascade:tosysobj:noStages'));
    end
    Hs = dsp.FilterCascade(sysobj(this.Stage(1)));
    for idx = 2:this.nstages
        Hs.addStage(sysobj(this.Stage(idx)));
    end
end

% --- Main helper functions ---

function setProperties(dfiltObj, sysObj)

% Parse dfiltObj and get information on branches, delays and gains
info = getInfoOnAllpassbranches(dfiltObj);

% Coefficients and delays in Branch #1. Branch1 #1 can be a pure delay
if(~info(1).IsDelay)
    % Branch not delay - assume info.Branch1.Filter is a
    % dfilt.cascade(wdf)allpass
    % Get coefficients for source dfilt in Branch1 as cell array
    C = getCoefficientsFromSourceDfilt(info(1));
    % If branch also includes delay, add further allpass stage to reflect
    % that. Cdelay is either empty or a cell array
    Cdelay = coefficientValuesForDelay(info(1));
    C = [C, Cdelay];
    % Set constructed Coefficients cell array in target sysObj, Branch #1
    setRelevantCoefficientsInTargetSysobj(C, sysObj, 1);
else
    % If Branch1 is simple delay, set PureDelayBranch to true and set the
    % Delay property to match the DelayLength found
    sysObj.PureDelayBranch = true;
    sysObj.Delay = info(1).DelayLength;
end

% Coefficients and delays for Branch #2 - this must contain a nontrivial
% cascaded allpass filter
% Get coefficients for source dfilt in Branch2 as cell array
C = getCoefficientsFromSourceDfilt(info(2));
% If branch also includes delay, add further allpass stage to reflect
% that. Cdelay is either empty or a cell array
Cdelay = coefficientValuesForDelay(info(2));
C = [C; Cdelay];

% Set relevant coefficients in System object (Branch #2), using available
% coefficients cell array
setRelevantCoefficientsInTargetSysobj(C, sysObj, 2);

% Set extra branch gains in target System object
setBranchGains(info, sysObj);

function warnIfStatesNontrivial(dfiltObj)
% This does not do anything, except throwing a warning if source dfilt
% object has populated internal states - Those will not be copied across to
% target System object

% To start, need to get to an actual nontrivial firlter object in the
% structure. Choose to go for the first reachable cascade(wdf)allpass
branchinfo = getEmptyInfo;

branchinfo = findFirstAllpassBranch(dfiltObj, branchinfo);

firstAllpassBranch = extractDfiltAtPath(dfiltObj, branchinfo.Path);

% Check if warning needed and throw if appropriate
if firstAllpassBranch.PersistentMemory
    % Get internal states
    DS = getinitialconditions(firstAllpassBranch);
    nontrivialStates = ~isempty(DS) && ~all(all(DS == 0));
    if(nontrivialStates)
        % The conversion to System object ignores
        % non-zero internal states as dsp.CoupledAllpassFilter does not
        % support non-zero initial conditions
        warning(message('signal:dfilt:cascade:tosysobj:nostatesupport'))
    end
end

% --- Secondary helper functions ---

function emptyinfo = getEmptyInfo

emptyinfo = struct(...
    'Path', [], ...
    'Filter', [], ...
    'IsDelay', [], ...
    'DelayLength', [], ...
    'Multiplier', []);

function info = getInfoOnAllpassbranches(idfilt)
% Parse idfilt and return information on the two parallel
% allpass branches in the structure, including
% - dfilt.cascade(wdf)allpass objects for each branch, where applicable
% - Whether or not each branch is a pure delay block and if yes, the
%   assaciated latency in sampling periods
% - The amound of additional delay to the allpass filter in the branch, if
%   any
% This assumes that
% - The overall structure of idfilt includes a dfilt.parallel
% - That dfilt.parallel has 2 stages
% - At least one of those stages includes (e.g. either directly or within a
%   further dfilt.cascade) a dfilt.cascade(wdf)allpass or a dfilt.dffir
% - If only one of the parallel stages includes a nontrivial dfilt object, 
%   then the other includes a dfilt.delay

% Find structural path of up to two cascade(wdf)allpass branches
info = getPathsOfAllpassBranches(idfilt);

% If two cascade(wdf)allpass branches are found, return collected info
if(isempty(info(2).Path))
    % Empty paths mean one allpass not found. Assume
    % - One cascade(wdf)allpass is always found
    % - If only one is found, then a delay is in other branch
    % If only one is found, keep it in Branch2 (make sure only Branch1 can
    % be set to pure delay
    tmp = info(2);
    info(2) = info(1);
    info(1) = tmp;

    % Get delay in Branch1 as parallel with cascade(wdf)allpas in Branch2
    delay1info = getDelayInParallelWithPath(idfilt, info(2));
    info(1).Path = delay1info.Path;
    info(1).IsDelay = true;
else
    % Get delay in Branch1 a sum of delaying objects in series with the
    % cascade(wdf)allpass filter identified
    delay1info = getDelayInSeriesWithPath(idfilt, info(1));
end
info(1).DelayLength = delay1info.DelayLength;

% Get delay in series with cascade(wdf)allpas in Branch2
delay2info = getDelayInSeriesWithPath(idfilt, info(2));
if(delay2info.IsDelay)
    info(2).DelayLength = delay2info.DelayLength;
end

% Get gains in series with paths
gain1info = getGainInSeriesWithPath(idfilt, info(1));
info(1).Multiplier = gain1info.Multiplier;

gain2info = getGainInSeriesWithPath(idfilt, info(2));
info(2).Multiplier = gain2info.Multiplier;

function info = getPathsOfAllpassBranches(idfilt)
% Returns [1x2] cell array with paths to up to two non-trivial (of type
% dsp.cascadeallpass, dsp.cascadewdfallpass or dfilt.dffir) leaf filters in
% overall filter structure

info = repmat(getEmptyInfo, 1, 2);

info(1) = findFirstAllpassBranch(idfilt, info(1));

tmp = copy(idfilt);

tmp = removeCascadeallpassFromPath(tmp, info(1).Path);

info(2) = findFirstAllpassBranch(tmp, info(2));

function oinfo = findFirstAllpassBranch(idfilt, iinfo)
% Returns a structural path pointing to the first cascade(wdf)allpass found
% in the structure hierarchy, as a vector of nested Stage numbers 
oinfo = iinfo;
ipath = iinfo.Path;
switch(class(idfilt))
    case {'dfilt.cascade', 'dfilt.parallel'}
        % Try Stage(1)
        oinfo.Path = [iinfo.Path, 1];
        oinfo.Filter = [];
        oinfo = findFirstAllpassBranch(idfilt.Stage(1), oinfo);
        if(isempty(oinfo.Path))
            % Try Stage(2)
            oinfo.Path = [iinfo.Path, 2];
            oinfo = findFirstAllpassBranch(idfilt.Stage(2), oinfo);
        end
    % case {'dfilt.cascadeallpass', 'dfilt.cascadewdfallpass', 'dfilt.dffir'}
    case {'dfilt.cascadeallpass', 'dfilt.cascadewdfallpass'}
        oinfo.Path = ipath;
        oinfo.Filter = idfilt;
        oinfo.IsDelay = false;
    otherwise % including dfilt.scalar, dfilt.dffir, dfilt.delay
        oinfo.Path = [];
        oinfo.Filter = [];
        
end

function odfilt = removeCascadeallpassFromPath(idfilt, path)
% Places a dfilt.scalar(0) at the location pointed to by the path vector

if(isempty(path))
    odfilt = idfilt;
    return
end

f = idfilt;
for k = 1:length(path)-1
    f = f.Stage(path(k));
end
f.Stage(path(end)) = dfilt.scalar(0);
odfilt = idfilt;

function sergaininfo = getGainInSeriesWithPath(idfilt, mainbranchinfo)

sergaininfo = getEmptyInfo;

path = mainbranchinfo.Path;

branchDfilt = extractParallelBranchContainingPath(idfilt,path);

gain = getTotalGainInCascade(branchDfilt);
sergaininfo.Multiplier = gain;

function serdelayinfo = getDelayInSeriesWithPath(idfilt, mainbranchinfo)

serdelayinfo = getEmptyInfo;

path = mainbranchinfo.Path;

branchDfilt = extractParallelBranchContainingPath(idfilt, path);

delaytaps = getTotalDelayInCascade(branchDfilt);

serdelayinfo.Path = path;
serdelayinfo.IsDelay = true;
serdelayinfo.DelayLength = delaytaps;

function branchcandidate = extractParallelBranchContainingPath(...
    dfiltObj,path)

branchcandidate = extractDfiltAtPath(dfiltObj, path);
for kh = length(path)-1:-1:1
    tmpdfilt = extractDfiltAtPath(dfiltObj, path(1:kh));
    if(isa(tmpdfilt, 'dfilt.parallel'))
        % We hit the top of where we want to search
        break
    end
    branchcandidate = tmpdfilt;
end

function pardelayinfo = getDelayInParallelWithPath(idfilt, mainbranchinfo)

pardelayinfo = getEmptyInfo;

path = mainbranchinfo.Path;

% Find parallel branch by walking upwards from given hierarchy path
for kh = length(path)-1:-1:1
    pdfilt = extractDfiltAtPath(idfilt, path(1:kh));
    if(isa(pdfilt, 'dfilt.parallel'))
        break
    end
end
if(isa(pdfilt, 'dfilt.parallel'))
    curBranch = path(kh + 1);
    parBranch = mod(curBranch,2)+1;
    % Here is the top-level parallel branch to the one passed as input
    parpath = [path(1:kh), parBranch];
    
    parbranch_dfilt = extractDfiltAtPath(idfilt, parpath);
    
    delaytaps = getTotalDelayInCascade(parbranch_dfilt);

    pardelayinfo.Path = parpath;
    % pardelayinfo.Filter = parbranch_dfilt;
    pardelayinfo.IsDelay = true;
    pardelayinfo.DelayLength = delaytaps;
end

function odfilt = extractDfiltAtPath(idfilt, path)
if(~isempty(path))
    odfilt = idfilt.Stage(path(1));
    for k = 2:length(path)
        odfilt = odfilt.Stage(path(k));
    end
else
    %     odfilt = [];
    odfilt = idfilt;
end

function delaytaps = getTotalDelayInCascade(varargin)
% getTotalDelayInCascade(indfilt, indelay)
assert(nargin == 1 || nargin == 2)
indfilt = varargin{1};

if(nargin == 2)
    indelay = varargin{2};
else
    indelay = 0;
end

switch(class(indfilt))
    case 'dfilt.cascade',
        % Sum delay over all cascade stages
        numstages = nstages(indfilt);
        delaytaps = indelay;
        for k = 1:numstages
            delaytaps = getTotalDelayInCascade(indfilt.Stage(k), delaytaps);
        end
    case 'dfilt.parallel',
        % Assuming to only be dealing with a cascade of dfilt.stages
        assert(false)
    case 'dfilt.scalar',
        delaytaps = indelay;
    case 'dfilt.delay',
        delaytaps = indelay + indfilt.Latency;
    case 'dfilt.dffir',
        delaytaps = indelay + median(grpdelay(indfilt, 3));
    otherwise
        % Assume no delay to report for other filter classes
        delaytaps = indelay;
end

function gain = getTotalGainInCascade(varargin)
% getTotalGainInCascade(indfilt, ingain)
assert(nargin == 1 || nargin == 2)
indfilt = varargin{1};

if(nargin == 2)
    ingain = varargin{2};
else
    ingain = 1;
end

switch(class(indfilt))
    case 'dfilt.cascade',
        % Sum delay over all cascade stages
        numstages = nstages(indfilt);
        gain = ingain;
        for k = 1:numstages
            gain = getTotalGainInCascade(indfilt.Stage(k), gain);
        end
    case 'dfilt.parallel',
        % Assuming to only be dealing with a cascade of dfilt.stages
        assert(false)
    case 'dfilt.scalar',
        gain = ingain * indfilt.Gain;
    case 'dfilt.delay',
        gain = ingain;
    case 'dfilt.dffir',
        % Only report gain if this has Numerator of type [0 ... 0 gain]
        if(all(indfilt.Numerator(1:end-1) == 0))
            gain = ingain * indfilt.Numerator(end);
        else
            gain = ingain;
        end
    otherwise
        % Assume no additional gain to report
        gain = ingain;
end

function C = getCoefficientsFromSourceDfilt(sourceBranchInfo)

CascadeFilter = sourceBranchInfo.Filter;
% Number of filter sections in top branch
numsec = length(fieldnames(CascadeFilter.AllpassCoefficients));
C = cell(numsec, 1);
for k = 1:numsec
    C{k} = CascadeFilter.AllpassCoefficients.(['Section',num2str(k)]);
    if(isa(CascadeFilter, 'dfilt.cascadewdfallpass'))
        C{k} = dsp.AllpassFilter.poly2wdf(C{k});
    end
end

function setRelevantCoefficientsInTargetSysobj(...
    Coefficients, sysObj, targetBranchNumber)

switch(sysObj.Structure)
    case 'Minimum multiplier'
        sysObj.(['AllpassCoefficients', num2str(targetBranchNumber)]) = ...
            Coefficients;
    case 'Wave Digital Filter'
        sysObj.(['WDFCoefficients', num2str(targetBranchNumber)]) = ...
            Coefficients;
    otherwise
        % Only filter classes dfilt.cascadeallpass and
        % dfilt.cascadewdfallpass expected
        assert(false)
end

function DelayCoefficients = ...
    coefficientValuesForDelay(sourceBranchInfo)

if(sourceBranchInfo.DelayLength > 0)
    if(isa(sourceBranchInfo.Filter, 'dfilt.cascadeallpass'))
        % The task is easy in this case as the 'Minimum multiplier'
        % Structure in dsp.AllpassFilter and dsp.CoupledAllpassFilter does
        % support arbitrary long single sections. Simply use one section
        DelayCoefficients = {zeros(1, sourceBranchInfo.DelayLength)};

    elseif(isa(sourceBranchInfo.Filter, 'dfilt.cascadewdfallpass'))
        % Things are more complicated in this case - the 'Wave Digital 
        % Filter' Structure in dsp.AllpassFilter and 
        % dsp.CoupledAllpassFilter only supports sections of order 1, 2 and
        % 4. Use as many 2nd-order sections as necessary plus possibly one
        % first-order section. All inner coefficients to be zeros

        % First take care of teh 2nd-order sections
        nSecondOrderSections = floor(sourceBranchInfo.DelayLength/2);
        DelayCoefficients = repmat({[0, 0]}, nSecondOrderSections, 1);
        % Add last 1st order section if necessary
        if(sourceBranchInfo.DelayLength - 2*nSecondOrderSections == 1)
            DelayCoefficients = [DelayCoefficients; {0}];
        end
    else
        % Only calsses dfilt.cascadeallpas and dfilt.cascadewdfallpas
        % expected as inputs
        assert(false)
    end
else
    DelayCoefficients = [];
end

function setBranchGains(info, sysObj)

switch(info(1).Multiplier)
    case 1,
        sysObj.Gain1 = '1';
    case -1,
        sysObj.Gain1 = '-1';
    case 1i,
        sysObj.Gain1 = '0+1i';
    case -1i,
        sysObj.Gain1 = '0-1i';
    otherwise
        error(message('signal:dfilt:cascade:tosysobj:badbranchgain'));
end
switch(info(2).Multiplier)
    case 1,
        sysObj.Gain2 = '1';
    case -1,
        sysObj.Gain2 = '-1';
    case 1i,
        sysObj.Gain2 = '0+1i';
    case -1i,
        sysObj.Gain2 = '0-1i';
    otherwise
        error(message('signal:dfilt:cascade:tosysobj:badbranchgain'));
end
