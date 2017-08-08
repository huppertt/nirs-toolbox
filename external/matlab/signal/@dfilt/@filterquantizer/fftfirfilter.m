function [y,z] = fftfirfilter(q,Hd,bfft,x,z)
%FFTFIRFILTER Filter this section.
%   [Y,Zf] = FFTFIRFILTER(q,Hd,bfft,X,ZI) filters this section.  This function is only
%   intended to be called from FFTFIRFILTER/SECFILTER.
%
%   See also DFILT.   
  
%   Author: R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.

x = quantizeinput(q,x);

L = Hd.BlockLength;
[Lx,nchannels] = size(x);
M = nstates(Hd);
Nfft = L+M;
y = zeros(Lx,nchannels); % Preallocate


if rem(Lx,L),    
    error(message('signal:dfilt:filterquantizer:fftfirfilter:InternalError'));
end

for n = 1:L:Lx,
    ytemp = ifft(bfft.*fft(x(n:n+L-1,:),Nfft));
    ytemp(1:M,:) = ytemp(1:M,:) + z;
    y(n:n+L-1,:) = ytemp(1:L,:);
    z = ytemp(end-M+1:end,:);
end


% We have removed scale values for now. 
% if Hd.ScaleValues~=1 && round(log2(Nfft))~=log2(Nfft),
%     warning('signal:dfilt:filterquantizer:fftfirfilter:scalevalues', 'Scale values ignored when the length of numerator + blocklength -1 is not a power of 2.');
% else
%     y = y*prod(Hd.ScaleValues);
% end

    
