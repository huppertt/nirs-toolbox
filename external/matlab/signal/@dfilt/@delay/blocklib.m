function [lib, srcblk,hasInputProcessing,hasRateOptions] = blocklib(~,~)
%BLOCKPARAMS Returns the library and source block for BLOCKPARAMS

% Copyright 2006-2012 The MathWorks, Inc.

b = isspblksinstalled;
if b,
    lib = 'dspobslib';
    srcblk = 'Delay';
else
    lib = 'simulink';
    srcblk = 'Discrete/Integer Delay';
end

hasInputProcessing = true;
hasRateOptions = false;
