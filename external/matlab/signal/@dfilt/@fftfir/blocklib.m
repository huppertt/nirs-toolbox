function [lib, srcblk,hasInputProcessing,hasRateOptions] = blocklib(~,~)
%BLOCKPARAMS Returns the library and source block for BLOCKPARAMS

% Copyright 1988-2012 The MathWorks, Inc.

% Library, block
lib = 'dsparch4';
srcblk = 'Overlap-Add FFT Filter';

hasInputProcessing = false;
hasRateOptions = false;
