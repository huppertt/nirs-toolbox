function pv = blockparams(this, mapstates, varargin) %#ok<INUSD>
%BLOCKPARAMS   Return the block parameters.

%   Copyright 1988-2012 The MathWorks, Inc.

% MAPSTATES is ignored because a gain doesn't have states.

pv  = scalarblockparams(this.filterquantizer);

pv.Gain = num2str(get(reffilter(this), 'Gain'));
