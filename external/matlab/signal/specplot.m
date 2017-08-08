function specplot(P,Fs)
%SPECPLOT Plot the output of the SPECTRUM function.
%   SPECPLOT will be removed in a future release of MATLAB.  Use one of the
%   following functions instead:
%      PERIODOGRAM
%      PWELCH
%      PBURG
%      PCOV
%      PMCOV
%      PYULEAR
%      PMTM
%      PEIG
%      PMUSIC
%      TFESTIMATE
%      MSCOHERE

%   SPECPLOT(P,Fs), uses P, the output of SPECTRUM, and Fs, the
%   sample frequency, to successively plot:
%
%   	Pxx - X Power Spectral Density & confidence.
%   	Pyy - Y Power Spectral Density & confidence.
%   	abs(Txy) - Transfer Function Magnitude.
%   	angle(Txy) - Transfer Function Phase.
%   	Cxy - Coherence Function.
%
%   The 95% confidence intervals are displayed on the power
%   spectral density curves.
%
%   SPECPLOT(P) uses normalized frequency, Fs = 2, so that 1.0 on
%   the frequency axis is half the sample rate (the Nyquist 
%   frequency).

%   Author(s): J.N. Little, 7-9-86
%   	   J.N. Little, 11-14-91, revised
%   Copyright 1988-2004 The MathWorks, Inc.

warning(message('signal:specplot:obsoleteFunction'));

[n,m] = size(P);
if nargin < 2
   Fs = 2;
end
f = (1:n-1)/n*Fs/2;

if m == 2
   c = [P(2:n,1)+P(2:n,2) P(2:n,1)-P(2:n,2)];
  else
   c = [P(2:n,1)+P(2:n,6) P(2:n,1)-P(2:n,6)];
end
c = c .* (c > 0);

newplot;
semilogy(f,P(2:n,1),f,c(:,1),'--',f,c(:,2),'--'), ...
title('Pxx - X Power Spectral Density'), ...
xlabel('Frequency')
if m == 2
   return
end

pause

newplot;
c = [P(2:n,2)+P(2:n,7) P(2:n,2)-P(2:n,7)];
c = c .* (c > 0);
semilogy(f,P(2:n,2),f,c(:,1),'--',f,c(:,2),'--'), ...
 title('Pyy - Y Power Spectral Density'), ...
 xlabel('Frequency'), pause

newplot;
semilogy(f,abs(P(2:n,4))), ...
 title('Txy - Transfer function magnitude'), ...
 xlabel('Frequency'), pause

newplot;
plot(f,180/pi*angle(P(2:n,4))), ...
 title('Txy - Transfer function phase'), ...
 xlabel('Frequency'), ...
 ylabel('Degrees'), pause

newplot;
plot(f,abs(P(2:n,5))), ...
 title('Cxy - Coherence'), ...
 xlabel('Frequency'), pause

