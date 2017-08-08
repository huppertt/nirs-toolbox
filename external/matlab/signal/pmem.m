function varargout = pmem( xR, p, nfft, Fs, flag )
%PMEM   Power Spectrum estimate via MEM (Maximum Entropy Method).
%   PMEM has been replaced by PYULEAR.  PMEM still works but may be
%   removed in the future.  Type help PYULEAR for details.
%
%   See also PBURG, PCOV, PMCOV, PMTM, PMUSIC, PEIG, PWELCH, PERIODOGRAM.

%   Ref: S.L. Marple, DIGITAL SPECTRAL ANALYSIS WITH APPLICATIONS,
%              Prentice-Hall, 1987, Chapter 7

%   Author(s): J. McClellan, 9-17-95
%   Copyright 1988-2012 The MathWorks, Inc.
%     

narginchk(2,5);

[Mx,Nx] = size(xR);
xIsMatrix = Mx>1 & Nx>1;

if nargin < 1
   error(message('signal:pmem:Nargchk'))
elseif  isempty(p)
   error(message('signal:pmem:Empty'))
end
if issparse(xR)
   error(message('signal:pmem:Sparse'))
end
if nargin < 5,   flag = [];   end
if nargin < 4,   Fs = [];     end  
if nargin < 3,   nfft = [];   end
if ischar(nfft)
   tflag = nfft; nfft = Fs; Fs = flag; flag = tflag;
elseif ischar(Fs)
   tflag = Fs;  Fs = flag; flag = tflag;
end

if isempty(nfft),    nfft = 256;  end
if isempty(Fs),      Fs = 2;      end
corr_flag = 0;     %<---- not expecting correlation
if ~isempty(flag)
   flag = upper(flag);
   if  (~isempty( strfind(flag,'CORR') )),  corr_flag = 1;  end
end

if( corr_flag )   %-- might be correlation matrix
   if Mx~=Nx
      error(message('signal:pmem:MatrixMustBeSquare'))
   elseif  norm(xR'-xR) > 100*eps
      error(message('signal:pmem:SignalErr'))
   end
end

if  ~xIsMatrix
   [RR,lags] = xcorr(xR,p,'biased'); %#ok
   xR = toeplitz( RR(p+1:2*p+1), RR(p+1:-1:1) );
else  %-- Matrix case
   if  p>= Mx
      error(message('signal:pmem:ColLengthMustBeGtOrder'))
   elseif  ~corr_flag
      xR = corrcoef( xR.' );
   end
end
a = [1;  xR(2:p+1,2:p+1)\(-xR(2:p+1,1)) ];
EE = abs( xR(1,1:p+1) * a );

Spec2 = abs( fft( a, nfft ) ) .^ 2;
Pxx = EE*ones(size(Spec2)) ./ Spec2;

%--- Select first half only, when input is real
if isreal(xR)   
    if rem(nfft,2),    % nfft odd
        select = (1:(nfft+1)/2)';
    else
        select = (1:nfft/2+1)';
    end
else
    select = (1:nfft)';
end

Pxx = Pxx(select);
ff = (select - 1)*Fs/nfft;

if nargout == 0
   newplot;
   plot(ff,10*log10(Pxx)), grid on
   xlabel('Frequency'), ylabel('Power Spectrum Magnitude (dB)');
   title('Maximum Entropy Spectral Estimate')
end

if nargout >= 1
    varargout{1} = Pxx;
end
if nargout >= 2
    varargout{2} = ff;
end
if nargout >= 3
    varargout{3} = a;
end
