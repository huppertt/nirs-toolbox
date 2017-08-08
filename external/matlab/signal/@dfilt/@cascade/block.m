function varargout = block(Hd, varargin)
%BLOCK Generate a DSP System Toolbox block equivalent to the filter object.
%   BLOCK(Hd, PARAMETER1, VALUE1, PARAMETER2, VALUE2, ...) allow you to
%   specify options in parameter/value pairs. Parameters can be:
%   'Destination': <'Current'>, 'New'
%   'Blockname': 'Filter' by default
%   'OverwriteBlock': 'on', <'off'>
%   'MapStates', 'on', <'off'>

%   Copyright 1988-2012 The MathWorks, Inc.

% Check if DSP System Toolbox is installed
[b, ~, ~, errObj] = isspblksinstalled;
if ~b
    error(errObj);
end

idx = find(strcmpi(varargin,'Link2Obj'));
if ~isempty(idx),
    link2obj = varargin{idx+1}; 
    if strcmpi(link2obj,'on'),
        error(message('signal:dfilt:cascade:block:noBlockLink'));
    end
end

% Check that all sections are supported
try
    hasInputProcessing = false(nstages(Hd),1);
    hasRateOptions = false(nstages(Hd),1);
    for i=1:nstages(Hd),
        [~, ~,hasInputProcessing(i),hasRateOptions(i)] = blocklib(Hd.Stage(i),'off');
        blockparams(Hd.Stage(i), 'off');
    end
catch ME
    error(message('signal:dfilt:cascade:block:NotSupported'));
end

% Parse inputs
[hTar, ~, ~, errObj]= uddpvparse('dspfwiztargets.blocktarget', varargin{:});
if ~isempty(errObj), error(errObj); end

% Create model
pos = createmodel(hTar);

% Creation of a subsystem
sys = hTar.system;
sysname = hTar.blockname;

if strcmpi(hTar.OverwriteBlock, 'on') %
    currentblk = find_system(sys, 'SearchDepth', 1,'LookUnderMasks', 'all', 'Name', sysname);
    if ~isempty(currentblk{1})
        delete_block(currentblk{1}); % Delete Filter block if present in the Destination
    end
end

xoffset = [100 0 100 0];
h = add_block('built-in/subsystem', sys, 'Tag', 'BlockMethodSubSystem');
if isempty(pos), pos = [65 40 140 80]; end
set_param(h,'Position',pos);

% Inport block
add_block('built-in/Inport', [sys '/In'], 'Position', [105 52 135 68]);
set_param(0, 'CurrentSystem', sys);
srcblk = 'In';

% Sections
mapstates = 'off';
idx = find(strcmpi(varargin,'MapStates'));
if ~isempty(idx), mapstates = varargin{idx+1}; end
        
% Map Coefficients to Ports 
% Determine coefficient names of filter in each stage and store the names
% in hTar.
try
    [hTar,doMapCoeffs2Ports] = parse_coeffstoexport(Hd,hTar);
catch ME
    throwAsCaller(ME);
end

pos = [65 40 140 80];

% Check stage by stage allowed input processing and rate option
% Default combination (columnsaschannels,enforcesinglerate) is supported by
% all filters.

% Check if input processing and rateoptions were specified
hTar.IsInputProcessingSpecified = false;
idx = find(strcmpi(varargin,'InputProcessing'), 1);
if ~isempty(idx)
  hTar.IsInputProcessingSpecified = true;
end
hTar.IsRateOptionSpecified = false;
idx = find(strcmpi(varargin,'RateOption'), 1);
if ~isempty(idx)
  hTar.IsRateOptionSpecified = true;
end

verifyinputprocrateopts(Hd,hTar, 'block');
    
for i=1:nstages(Hd),
    secname = ['Stage' sprintf('%d',i)];
    if doMapCoeffs2Ports
        seccoeffnames = hTar.CoeffNames.(sprintf('Stage%d',i));
        params = {'Blockname', secname, 'MapStates', mapstates,...
                'MapCoeffsToPorts','on','CoeffNames',seccoeffnames};        
    else
       params = {'Blockname', secname, 'MapStates', mapstates};
    end
    
    if hasInputProcessing(i)
      params = [params {'InputProcessing', hTar.InputProcessing}]; %#ok<*AGROW>
    end
    if hasRateOptions(i)
      params = [params {'RateOption', hTar.RateOption}];
    end
    
    block(Hd.Stage(i),params{:});
             
    set_param(0, 'CurrentSystem', sys);
    set_param([sys, '/', secname], 'Position', pos+i*xoffset);
    add_line(sys,[srcblk '/1'], [secname '/1'], 'autorouting', 'on');
    srcblk = secname;
end

% Outport block
outblk = add_block('built-in/Outport', [sys '/Out']);
set_param(outblk, 'Position', [65 52 95 68]+(nstages(Hd)+1)*xoffset)
add_line(sys,[srcblk '/1'], 'Out/1', 'autorouting', 'on');

if nargout,
    varargout = {h};
end

% Refresh connections
oldpos = get_param(sys, 'Position');
set_param(sys, 'Position', oldpos + [0 -5 0 -5]);
set_param(sys, 'Position', oldpos);

% Open system
slindex = strfind(sys,'/');

% When the model is not reused and when the path points to a system or a
% non-mask block, open the system
if isMaskOff(hTar)
    % Open system
    open_system(sys(1:slindex(end)-1));
end

if nargout,
    varargout = {h};
end

% [EOF]
