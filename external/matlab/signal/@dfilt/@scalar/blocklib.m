function [lib, srcblk,hasInputProcessing,hasRateOptions] = blocklib(~,~)
%BLOCKPARAMS Returns the library and source block for BLOCKPARAMS

% Copyright 2006-2012 The MathWorks, Inc.

lib = 'built-in';
srcblk = 'Gain';

hasInputProcessing = false;
hasRateOptions = false;
