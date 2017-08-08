function [lib, srcblk,hasInputProcessing,hasRateOptions] = blocklib(Hd,~)
%BLOCKPARAMS Returns the library and source block for BLOCKPARAMS

%   Copyright 2008-2012 The MathWorks, Inc.

% Library, block

lib = 'dsparch4';

checksv(Hd)
srcblk = 'Biquad Filter';

hasInputProcessing = true;
hasRateOptions = false;
