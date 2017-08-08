function [N, Wn, bta, filtype] = kaiserord(fcuts, mags, devs, fsamp, cellflag) %#ok
%KAISERORD FIR order estimator (lowpass, highpass, bandpass, multiband).
%   [N,Wn,BTA,FILTYPE] = KAISERORD(F,A,DEV,Fs) is the approximate order N, 
%   normalized frequency band edges Wn, Kaiser window beta parameter BTA
%   and filter type FILTYPE to be used by the FIR1 function:
%      B = FIR1(N, Wn, FILTYPE, kaiser( N+1,BTA ), 'noscale' )
%
%   The resulting filter will approximately meet the specifications given
%   by the input parameters F, A, and DEV.
%
%   F is a vector of band edge frequencies in Hz, in ascending order
%   between 0 and half the sampling frequency Fs.  A is a vector of 0s and
%   1s specifying the desired function's amplitude on the bands defined by
%   F. The length of F is twice the length of A, minus 2 (it must therefore
%   be even).  The first frequency band is assumed to start at zero, and
%   the last one always ends at Fs/2.
%
%   DEV is a vector of maximum deviations or ripples (in linear units)
%   allowable for each band. The smallest deviation specified (MIN(DEV)) is
%   used for both the passband and the stopband.
%
%   Fs is the sampling frequency (which defaults to 2 if you leave it off).
%
%   C = KAISERORD(F,A,DEV,Fs,'cell') is a cell-array whose elements are the
%   parameters to FIR1.
%
%   EXAMPLE
%      Design a lowpass filter with a passband edge of 1500Hz, a 
%      stopband edge of 2000Hz, passband ripple of 0.01, stopband ripple 
%      of 0.1, and a sampling frequency of 8000Hz:
%
%      [n,Wn,bta,filtype] = kaiserord( [1500 2000], [1 0], [0.01 0.1], 8000 );
%      b = fir1(n, Wn, filtype, kaiser(n+1,bta), 'noscale');
%   
%      This is equivalent to
%      c = kaiserord( [1500 2000], [1 0], [0.01 0.1], 8000, 'cell' );
%      b = fir1(c{:});
%
%   CAUTION 1: The order N is just an estimate. If the filter does not
%   meet the original specifications, a higher order such as N+1, N+2, etc. 
%   will; if the filter exceeds the specs, a slightly lower order one may work.
%   CAUTION 2: Results are inaccurate if cutoff frequencies are near zero
%   frequency or the Nyquist frequency; or if the devs are large (10%).
%
%   See also FIR1, KAISER, FIRPMORD.

%   Author(s): J. H. McClellan, 10-28-91
%   Copyright 1988-2013 The MathWorks, Inc.

%   References:
%  [1] J.F. Kaiser, ``Nonrecursive Digital Filter Design Using
%       the I_o-sinh Window Function,'' Proc. 1974 IEEE
%       Symp. Circuits and Syst., April 1974, pp. 20--23.
%  [2] IEEE, Digital Signal Processing II, IEEE Press, New York:
%      John Wiley & Sons, 1975, pp. 123--126.

narginchk(3,5)
if nargin < 4 || isempty(fsamp),
    fsamp = 2;
end

% Cast to enforce precision rules
fcuts = signal.internal.sigcasttofloat(fcuts,'double','kaiserord','F',...
  'allownumeric');
mags = signal.internal.sigcasttofloat(mags,'double','kaiserord','A',...
  'allownumeric');
devs = signal.internal.sigcasttofloat(devs,'double','kaiserord','DEV',...
  'allownumeric');
fsamp = signal.internal.sigcasttofloat(fsamp,'double','kaiserord','Fs',...
  'allownumeric');

if max(fcuts) >= fsamp/2,
    error(message('signal:kaiserord:InvalidRange'));
end


fcuts = fcuts/fsamp;       %  NORMALIZE to sampling frequency

% Turn vectors into column vectors
fcuts = fcuts(:);
mags = mags(:);
devs = devs(:);

mf = size(fcuts,1);
nbands = size(mags,1);

if size(mags,1) ~= size(devs,1)
    error(message('signal:kaiserord:InvalidDimensionsADEV', 'A', 'DEV'));
end
if( min(abs(mags)) )
   error(message('signal:kaiserord:SignalErrStopbands'));
end
dmags = abs(diff(mags));
if( any(dmags~=dmags(1)) )
    error(message('signal:kaiserord:SignalErrPassbands'));
end
if( any(diff(fcuts)<0) )
    error(message('signal:kaiserord:InvalidFreqVec'));
end


if mf ~= 2*(nbands-1)
    error(message('signal:kaiserord:InvalidDimensionsLengthF', 'F', '2*length(A)-2'));
end

zz = mags==0;             % find stopbands
devs = devs./(zz+mags);   % divide delta by mag to get relative deviation

% Determine the smallest width transition zone
% Separate the passband and stopband edges
%
f1 = fcuts(1:2:(mf-1));
f2 = fcuts(2:2:mf);
[df,n] = min(f2-f1); %#ok

%=== LOWPASS case: Use formula (ref: Herrmann, Rabiner, Chan)
%
if( nbands==2 )
     [L,bta] = kaislpord( f1(n), f2(n), devs(1), devs(2));

%=== BANDPASS case:
%    - try different lowpasses and take the WORST one that
%        goes through the BP specs; try both transition widths
%    - will also do the bandreject case
%    - does the multi-band case, one bandpass at a time.
%    
else
  L = 0;  bta = 0;
  for i=2:nbands-1,
    [L1,bta1] = kaislpord( f1(i-1), f2(i-1), devs(i),   devs(i-1) );
    [L2,bta2] = kaislpord( f1(i),   f2(i),   devs(i),   devs(i+1) );
    if( L1>L )
        bta = bta1;  L = L1;   end
    if( L2>L )
        bta = bta2;  L = L2;   end
  end
end

N = ceil( L ) - 1;   % need order, not length, for Filter design

%=== Make the MATLAB compatible specs for FIR1
%
Wn = 2*(f1+f2)/2;    %-- use mid-frequency; multiply by 2 for MATLAB
filtype = 'low';
if( nbands==2 && mags(1)==0 )
  filtype='high';
elseif( nbands==3 && mags(2)==0 )
  filtype='stop';
elseif( nbands>=3 && mags(1)==0 )  
  filtype='DC-0';                    
elseif( nbands>=3 && mags(1)==1 ) 
  filtype='DC-1';                   
end

% If order is odd, and gain is not zero at nyquist, increase the order by one.
if rem(N,2) && mags(end)~=0,
    N = N + 1;
end

if nargout == 1 && nargin == 5
  N = {N, Wn, filtype, kaiser(N+1,bta), 'noscale'};
end

%%%% ---- end of kaiserord

function [L,bta] = kaislpord(freq1, freq2, delta1, delta2 )
%KAISLPORD FIR lowpass filter Length estimator
%
%   [L,bta] = kaislpord(freq1, freq2, dev1, dev2)
%
%   input:
%     freq1: passband cutoff freq (NORMALIZED)
%     freq2: stopband cutoff freq (NORMALIZED)
%      dev1: passband ripple (DESIRED)
%      dev2: stopband attenuation (not in dB)
%
%   outputs:
%      L = filter Length (# of samples)   **NOT the order N, which is N = L-1
%   bta =  parameter for the Kaiser window
%
%   NOTE: Will also work for highpass filters (i.e., f1 > f2)
% 	      Will not work well if transition zone is near f = 0, or
%         near f = fs/2

%
% 	Author(s): J. H. McClellan, 8-28-95
	
%   References:
%     [1] Rabiner & Gold, Theory and Applications of DSP, pp. 156-7.     

delta = min( [delta1,delta2] );
atten = -20*log10( delta );
D = (atten - 7.95)/(2*pi*2.285);   %--- 7.95 was in Kaiser's original paper
%
df = abs(freq2 - freq1);
%
L = D/df + 1;
%
bta = signal.internal.kaiserBeta(atten);
