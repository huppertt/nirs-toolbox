function [Pxx, F, Pxxc] = psdcenterdc(Pxx, F, Pxxc, psdoptions)
%PSDCENTERDC  Center power and frequency for a given psdoptions structure
%   [PXX, F] = PSDCENTERDC(PXX,F,PSDOPTIONS) centers the power spectrum and
%   frequency vector based upon the Fs, nfft, and range fields in
%   PSDOPTIONS
 
%   Copyright 2014 The MathWorks, Inc.

nFreq = numel(F);
if nFreq == 0
  return
end

iseven = psdoptions.nfft/2 == round(psdoptions.nfft/2);
isonesided = strcmpi(psdoptions.range,'onesided');

if isonesided
  % Undo any x2 scaling of frequencies and confidence interval estimates
  if iseven 
    if ~isfield(psdoptions,'eigenvals')
      % divide all powers by 2 except nyquist and DC
      Pxx(2:end-1,:) = Pxx(2:end-1,:)/2;
      if ~isempty(Pxxc)
        Pxxc(2:end-1,:) = Pxxc(2:end-1,:)/2;
      end
    end
    idx = [nFreq-1:-1:2 1:nFreq];
  else
    if ~isfield(psdoptions,'eigenvals')
      % divide all powers by 2 except DC
      Pxx(2:end) = Pxx(2:end,:)/2;
      if ~isempty(Pxxc)
        Pxxc(2:end,:) = Pxxc(2:end,:)/2;
      end
    end
    idx = [nFreq:-1:2 1:nFreq];
  end
else
  if iseven
    idx = [nFreq/2+2:nFreq 1:nFreq/2+1];
  else
    idx = [(nFreq+1)/2+1:nFreq 1:(nFreq+1)/2];
  end
end

Pxx = Pxx(idx,:);
F = F(idx);
if ~isempty(Pxxc)
  Pxxc = Pxxc(idx,:);
end

Fs = psdoptions.Fs;
if isempty(Fs)
  % normalize to 2*pi when default specified.
  Fs = 2*pi;
end

if isonesided
  F(1:(end-nFreq)) = -F(1:(end-nFreq));
elseif iseven
  F(1:nFreq/2-1) = F(1:nFreq/2-1) - Fs;
else
  F(1:(nFreq-1)/2) = F(1:(nFreq-1)/2) - Fs;
end


