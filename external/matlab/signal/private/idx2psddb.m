function [psd, freq] = idx2psddb(Pxx, F, idx)
%IDX2PSDDB report back 10*log10(Pxx) and frequency at specified index
%   return NaN when index is empty
%   
%   This function is for internal use only. It may be removed in the future.
   
%   Copyright 2013 The MathWorks, Inc.

if ~isempty(idx)
  psd = 10*log10(Pxx(idx));
  freq = F(idx);
else
  psd = NaN;
  freq = NaN;
end
