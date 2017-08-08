function [bz,az] = impinvar(b,a,Fs,tol)
%IMPINVAR Impulse invariance method for analog to digital filter conversion.
%   [BZ,AZ] = IMPINVAR(B,A,Fs) creates a digital filter with numerator
%   and denominator coefficients BZ and AZ respectively whose impulse
%   response is equal to the impulse response of the analog filter with
%   coefficients B and A sampled at a frequency of Fs Hertz.  The B and A
%   coefficients will be scaled by 1/Fs.
%
%   If you don't specify Fs, it defaults to 1 Hz.
%
%   [BZ,AZ] = IMPINVAR(B,A,Fs,TOL) uses the tolerance TOL for grouping
%   repeated poles together.  Default value is 0.001, i.e., 0.1%.
%
%   NOTE: the repeated pole case works, but is limited by the
%         ability of the function ROOTS to factor such polynomials.
%
%   % EXAMPLE: Illustrate the relationship between digital and analog
%   % frequency responses.
%   [b,a] = butter(4,0.3,'s');
%   [bz,az] = impinvar(b,a,10);
%   [Ha,Wa] = freqs(b,a,512);  [Hz,Wz] = freqz(bz,az,512,10);
%   plot(Wa/(2*pi),20*log10(abs(Ha)),'LineWidth',2); hold on;
%   plot(Wz,20*log10(abs(Hz)),'r--');
%   xlabel('Frequency (Hz)'), ylabel('Magnitude (dB)');
%   title('Magnitude Response Comparison');
%   legend('Analog Filter','Digital Filter');
%
%   See also BILINEAR.

%   Added multiple pole capability, 3-Feb-96  J McClellan
%   Also, the filter gain is now multiplied by 1/Fs (per O&S, etc)

%   Author(s): J. McClellan, Georgia Tech, EE, DSP, 1990
%   Copyright 1988-2013 The MathWorks, Inc.

narginchk(2,4)
% Cast to enforce precision rules
if (any([signal.internal.sigcheckfloattype(b,'single','impinvar','B')...
    signal.internal.sigcheckfloattype(a,'single','impinvar','A')]))
    b = single(b);
    a = single(a);
end

if nargin<4,  tol = 1e-3; end
if nargin<3,  Fs = 1; end
% Cast to enforce precision rules
Fs = signal.internal.sigcasttofloat(Fs,'double','impinvar','Fs',...
  'allownumeric');
tol = signal.internal.sigcasttofloat(tol,'double','impinvar','TOL',...
  'allownumeric');

[M,N] = size(a);
if ~xor(M==1,N==1)
    error(message('signal:impinvar:MustBeVector', 'A'))
end
[M,N] = size(b);
if M > 1 && N > 1
    error(message('signal:impinvar:MustBeVector', 'B'))
end

b = b(:);
a = a(:);
a1 = a(1);
if a1 ~= 0
    % Ensure monotonicity of a
    a = a/a1;
end
kimp=[];
if length(b) > length(a)
    error(message('signal:impinvar:InvalidRange'))
elseif  (length(b)==length(a))
    % remove direct feed-through term, restore later
    kimp = b(1)/a(1);
    b = b - kimp*a;  b(1)=[];
end

%--- Achilles Heel is repeated roots, so I adapted code from residue
%---  and resi2 here.  Now I can group roots, and there is no division.
pt = roots(a).';
Npoles = length(pt);
[mm,ip] = mpoles(pt,tol);
pt = pt(ip);
starts = find(mm==1);
ends = [starts(2:length(starts))-1;Npoles];
for k = 1:length(starts)
    jkl = starts(k):ends(k);
    polemult(jkl) = mm(ends(k))*ones(size(jkl)); %#ok
    poleavg(jkl) = mean(pt(jkl))*ones(size(jkl)); %#ok
end
rez = zeros(Npoles,1);
kp = Npoles;
while kp>0
    pole = poleavg(kp);
    mulp = polemult(kp);
    num = b;
    den = poly( poleavg([1:kp-mm(kp),kp+1:Npoles]) );
    rez(kp) = polyval(num,pole) ./ polyval(den,pole);
    kp = kp-1;
    for k=1:mulp-1
        [num,den] = polyder(num/k,den);
        rez(kp) = polyval(num,pole) ./ polyval(den,pole);
        kp = kp-1;
    end
end

%----- Now solve for H(z) residues via impulse response matching
r = rez./gamma(mm);
p = poleavg;

az = poly(exp(p/Fs)).';
tn = (0:Npoles-1)'/Fs;
mm1 = mm(:).' - 1;
tt = tn(:,ones(1,Npoles)) .^ mm1(ones(size(tn)),:);
ee = exp(tn*p);
hh = ( tt.*ee ) * r;
bz = filter(az,1,hh);

if ~isempty(kimp)
    % restore direct feed-through term
    bz = kimp*az(:) + [bz(:);0];
end
bz = bz/Fs;

bz = bz(:).';   % make them row vectors
az = az(:).';

cmplx_flag = any(imag(b)) | any(imag(a));
if ~cmplx_flag
    if  norm(imag([bz az]))/norm([bz az]) > 1000*eps
        warning(message('signal:impinvar:NonRobust'));
    end
    bz = real(bz);
    az = real(az);
end
if a1~=0
    az = az*a1;
end
