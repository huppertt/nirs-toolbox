function s = blockparams(Hd, mapstates, varargin)
%BLOCKPARAMS Returns the parameters for BLOCK

% Copyright 1988-2012 The MathWorks, Inc.


% Parameters of the block
Nfft = Hd.BlockLength+Hd.ncoeffs-1;
B_Nfft = 2^nextpow2(Nfft);
if B_Nfft==Nfft,
    s.Nfft = num2str(Nfft);
else
    suggested_BlockLength = B_Nfft-Hd.ncoeffs+1;
    error(message('signal:dfilt:fftfir:blockparams:InvalidNFFT', num2str( suggested_BlockLength )));
end
s.h = mat2str(Hd.Numerator,18);

if ~isreal(Hd.Numerator),
    s.output_complexity = 'Complex';
end

% IC
if strcmpi(mapstates, 'on'),
    warning(message('signal:dfilt:fftfir:blockparams:mappingstates', srcblk));
end

