function pv = blockparams(this, mapstates, varargin)
%BLOCKPARAMS   Return the block parameters.

%   Copyright 1988-2012 The MathWorks, Inc.

b = isspblksinstalled;
if b,
    % Use DSP System Toolbox block
    pv.delay = sprintf('%d',this.Latency);
else
    pv.NumDelays = sprintf('%d',this.Latency);
end

% IC
if strcmpi(mapstates, 'on'),
    pv.IC = mat2str(getinitialconditions(this));
end

