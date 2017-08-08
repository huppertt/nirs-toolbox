function [N, ff, aa, wts] = firpmord(fcuts, mags, devs, fsamp, cellflag) %#ok
%FIRPMORD  Parks-McClellan optimal equiripple FIR order estimator.
%   [N,Fo,Ao,W] = FIRPMORD(F,A,DEV,Fs) finds the approximate order N, 
%   normalized frequency band edges Fo, frequency band magnitudes Ao and 
%   weights W to be used by the FIRPM function as follows:
%       B = FIRPM(N,Fo,Ao,W)
%   The resulting filter will approximately meet the specifications given
%   by the input parameters F, A, and DEV.  F is a vector of cutoff
%   frequencies in Hz, in ascending order between 0 and half the sampling
%   frequency Fs. If you do not specify Fs, it defaults to 2.  A is a
%   vector specifying the desired function's amplitude on the bands defined
%   by F. The length of F is twice the length of A, minus 2 (it must
%   therefore be even).  The first frequency band always starts at zero,
%   and the last always ends at Fs/2.  It is not necessary to add these
%   elements to the  F vector.  DEV is a vector of maximum deviations or
%   ripples (in linear units) allowable for each band.  DEV must have the
%   same length as A.
%
%   C = FIRPMORD(F,A,DEV,FS,'cell') is a cell-array whose elements are the 
%   parameters to FIRPM.
%
%   EXAMPLE
%      Design a lowpass filter with a passband-edge frequency of 1500Hz, a 
%      stopband-edge of 2000Hz, passband ripple of 0.01, stopband ripple 
%      of 0.1, and a sampling frequency of 8000Hz:
%
%         [n,fo,mo,w] = firpmord( [1500 2000], [1 0], [0.01 0.1], 8000 );
%         b = firpm(n,fo,mo,w);
%
%      This is equivalent to
%         c = firpmord( [1500 2000], [1 0], [0.01 0.1], 8000, 'cell');
%         b = firpm(c{:});
%
%   CAUTION 1: The order N is often underestimated. If the filter does not
%   meet the original specifications, a higher order such as N+1 or N+2 will.
%   CAUTION 2: Results are inaccurate if cutoff frequencies are near zero
%   frequency or the Nyquist frequency.
%
%   See also FIRPM, KAISERORD, DESIGNFILT.

%   Author(s): J. H. McClellan, 10-28-91
%   Copyright 1988-2013 The MathWorks, Inc.
%     

%   References:
%     [1] Rabiner & Gold, Theory and Applications of DSP, pp. 156-7.     

narginchk(3,5)
if nargin == 3,
    fsamp = 2;
end
% Cast to enforce precision rules
fcuts = signal.internal.sigcasttofloat(fcuts,'double','firpmord','F',...
  'allownumeric');
mags = signal.internal.sigcasttofloat(mags,'double','firpmord','A',...
  'allownumeric');
devs = signal.internal.sigcasttofloat(devs,'double','firpmord','DEV',...
  'allownumeric');
fsamp = signal.internal.sigcasttofloat(fsamp,'double','firpmord','Fs',...
  'allownumeric');

fcuts = fcuts/fsamp;       %  NORMALIZE to sampling frequency
if any(fcuts > 1/2),
    error(message('signal:firpmord:InvalidRange'));
end
if any(fcuts < 0),
    error(message('signal:firpmord:MustBePositive'));
end

% Turn vectors into column vectors
fcuts = fcuts(:);
mags = mags(:);
devs = devs(:);

mf = size(fcuts, 1);
mm = size(mags, 1);
nbands = mm;

if length(mags) ~= length(devs)
    error(message('signal:firpmord:MismatchedVectorLength', 'A', 'DEV'))
end

if mf ~= 2*(nbands-1)
    error(message('signal:firpmord:InvalidLength', 'F', '2*length(A)-2'))
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
     L = remlpord( f1(n), f2(n), devs(1), devs(2));

%=== BANDPASS case:
%    - try different lowpasses and take the WORST one that
%        goes through the BP specs; try both transition widths
%    - will also do the bandreject case
%    - does the multi-band case, one bandpass at a time.
%    
else
  L = 0;
  for i=2:nbands-1,
    L1 = remlpord( f1(i-1), f2(i-1), devs(i),   devs(i-1) );
    L2 = remlpord( f1(i),   f2(i),   devs(i),   devs(i+1) );
    L = max( [L; L1; L2] );
  end
end

N = ceil( L ) - 1;   % need order, not length, for firpm

%=== Make the MATLAB compatible specs for firpm
%
ff = [0;2*fcuts;1];
am(1:2:2*nbands-1) = mags;	
aa = [am';0] + [0;am'];
wts = ones(size(devs))*max(devs)./devs;

% If gain is not zero at nyquist, the order must be even.
% If the order is odd, we bump up the order by one.
if (aa(end) ~= 0) && rem(N,2),
    N = N+1;
end

if nargout == 1 && nargin == 5
   N = {N, ff, aa, wts};
end

%----------------------------------------------------------------------------
function L = remlpord(freq1, freq2, delta1, delta2 )
%REMLPORD FIR lowpass filter Length estimator
%
%   L = remlpord(freq1, freq2, dev1, dev2)
%
%   input:
%     freq1: passband cutoff freq (NORMALIZED)
%     freq2: stopband cutoff freq (NORMALIZED)
%     dev1:  passband ripple (DESIRED)
%     dev2:  stopband attenuation (not in dB)
%
%   output:
%     L = filter Length (# of samples)
%         NOT the order N, which is N = L-1
%
%   NOTE: Will also work for highpass filters (i.e., f1 > f2)
%         Will not work well if transition zone is near f = 0, or
%         near f = fs/2

% 	Author(s): J. H. McClellan, 10-28-91
%       Was Revision: 1.4, Date: 1994/01/25 17:59:46
	
%   References:
%     [1] Rabiner & Gold, Theory and Appications of DSP, pp. 156-7.     

AA = [-4.278e-01  -4.761e-01  0;
      -5.941e-01   7.114e-02  0;
      -2.660e-03   5.309e-03  0 ];
d1 = log10(delta1);
d2 = log10(delta2);
D = [1 d1 d1*d1] * AA * [1; d2; d2*d2];
%
bb = [11.01217; 0.51244];
fK =  [1.0  (d1-d2)] * bb;
%
df = abs(freq2 - freq1);
%
L = D/df - fK*df + 1;

% [EOF] - FIRPMORD.M
