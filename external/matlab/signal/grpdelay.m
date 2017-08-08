function [gd_out,w_out] = grpdelay(b,varargin)
%GRPDELAY Group delay of digital filter
%   [Gd,W] = GRPDELAY(B,A,N) returns length N vectors Gd and W containing
%   the group delay, and the frequencies (in radians) at which it is
%   evaluated, for the filter:
%
%               jw               -jw              -jmw
%        jw  B(e)    b(1) + b(2)e + .... + b(m+1)e
%     H(e) = ---- = ------------------------------------
%               jw               -jw              -jnw
%            A(e)    a(1) + a(2)e + .... + a(n+1)e
%
%   given numerator and denominator coefficients in vectors B and A.
%
%   [Gd,W] = GRPDELAY(SOS,N) computes the group delay of the filter
%   specified using the second order sections matrix SOS. SOS is a Kx6
%   matrix, where the number of sections, K, must be greater than or equal
%   to 2. Each row of SOS corresponds to the coefficients of a second order
%   filter. From the transfer function displayed above, the ith row of the
%   SOS matrix corresponds to [bi(1) bi(2) bi(3) ai(1) ai(2) ai(3)].
%
%   [Gd,W] = GRPDELAY(D,N) computes the group delay of the digital filter,
%   D. You design a digital filter, D, by calling the <a href="matlab:help designfilt">designfilt</a> function.
%
%   Group delay is defined as -d{angle(w)}/dw. The frequency response is
%   evaluated at N points equally spaced around the upper half of the unit
%   circle. If you don't specify N, it defaults to 512.
%
%   [Gd,W] = GRPDELAY(...,N,'whole') uses N points around the whole unit
%   circle.
%
%   [Gd,F] = GRPDELAY(...,N,Fs) and [Gd,F] = GRPDELAY(...,N,'whole',Fs)
%   return a frequency vector, F, in Hz when you specify the sample rate Fs
%   in Hz.
%
%   Gd = GRPDELAY(...,W) and Gd = GRPDELAY(..,F,Fs) return the group delay
%   evaluated at the points specified in frequency vectors W (in
%   radians/sample), or F (in Hz).
%
%   GRPDELAY(...) with no output arguments plots the group delay in the
%   current figure window.
%
%   % Example 1:
%   %   Design a lowpass FIR filter with normalized cut-off frequency at 
%   %   0.3 and determine its group delay.
%
%   b=fircls1(54,0.3,0.02,0.008);   
%   grpdelay(b)                     
%
%   % Example 2: 
%   %   Design a 5th order lowpass elliptic IIR filter and determine its  
%   %   group delay response.
%
%   [b,a] = ellip(5,0.5,20,0.4);    
%   grpdelay(b,a,512,'whole');      
%
%   % Example 3:
%   %   Design a Butterworth highpass IIR filter, represent its coefficients
%   %   using second order sections, and display its group delay response.
%
%   [z,p,k] = butter(6,0.7,'high');
%   SOS = zp2sos(z,p,k);    
%   grpdelay(SOS)                   
%
%   % Example 4:
%   %   Use the designfilt function to design a highpass IIR digital filter 
%   %   with order 8, passband frequency of 75 KHz, and a passband ripple 
%   %   of 0.2 dB. Sample rate is 200 KHz. Visualize the group delay response 
%   %   using 2048 frequency points.
%  
%   D = designfilt('highpassiir', 'FilterOrder', 8, ...
%            'PassbandFrequency', 75e3, 'PassbandRipple', 0.2,...
%            'SampleRate', 200e3);
%
%   grpdelay(D,2048)
%
%   See also FREQZ.

% Group delay algorithm notes:
%
% The Smith algorithm is employed for FIR filters
%  - unpublished algorithm from J. O. Smith, 9-May-1988
%  - faster than the Shpak algorithm, best suited to FIR filters
%
% The Shpak algorithm is employed for IIR filters
%  - it is more accurate than Smith algorithm for "tricky" IIR filters
%  - algorithm converts filter to second-order sections before computation
%  - D. Shpak, 17-Feb-2000
%  - Added special treatment of singularities based on a suggestion
%    from P. Kabal and extended to complex polys.  D. Shpak, 24-Nov-2001
%
%   Copyright 1988-2013 The MathWorks, Inc.

narginchk(1,5)
inputArgs = {};

isTF = true; % True if dealing with a transfer function

if all(size(b)>[1 1])
  % Input is a matrix, check if it is a valid SOS matrix
  if size(b,2) ~= 6
    error(message('signal:signalanalysisbase:invalidinputsosmatrix'));
  end
  isTF = false; % SOS instead of transfer function  
  
  % add a dummy input to varargin so that we can parse the rest of the
  % inputs with the same code
  varargin = [{1} varargin];  
end


% Define parameter defaults:
isWholeUnitCircle = 0;
isNormalizedFreq = 1;
Fs = 1;

vararginLength = length(varargin);
% Parse user arguments
if vararginLength == 0
    a = 1;  
    n = 512;
elseif vararginLength == 1
    a = varargin{1};    
    if all(size(a)>[1 1])
      error(message('signal:signalanalysisbase:inputnotsupported'));
    end        
    n = 512;
elseif vararginLength == 2
    a = varargin{1};
    n = varargin{2};
elseif vararginLength == 3
    a = varargin{1};
    n = varargin{2};
    dum = varargin{3};
    if ischar(dum),
        if strcmpi(dum, 'whole'),
            isWholeUnitCircle = 1;
            inputArgs = [inputArgs {'whole'}];
        end
    else
        isNormalizedFreq = 0;
        Fs = dum;
        inputArgs = [inputArgs {Fs}];
    end
elseif vararginLength == 4
    a = varargin{1};
    n = varargin{2};
    dum = varargin{3};
      
    if ~ischar(dum),
        error(message('signal:grpdelay:MustBeAString'));
    end    
    if strcmpi(dum, 'whole'),
        isWholeUnitCircle = 1;
        inputArgs = [inputArgs {'whole'}];
    end
    Fs = varargin{4};
    inputArgs = [inputArgs {Fs}];
    isNormalizedFreq = 0;
end

% Cast inputs to double. Check if any input is single so that we can cache
% the output to single precision. 
if ~all(cellfun(@isfloat,{b,a,n}))
  error(message('signal:signalanalysisbase:invalidDataType'));
end
isSingleArith = isa(b,'single') || isa(a,'single');

% Cast to enforce precision rules
n = double(n);
Fs = signal.internal.sigcasttofloat(Fs,'double','grpdelay','Fs',...
  'allownumeric');

if isTF
  % Choose group-delay computation algorithm:
  %
  % The Smith method is used for FIR filters
  %  - it is faster than the Shpak method and still accurate
  % The Shpak method is used for IIR filters
  %  - it is more accurate (converts to second-order sections)  
  
  % Determine if input is FIR
  if signalpolyutils('isfir',b,a)
    dlyFcn = @smithDly;
  else
    dlyFcn = @shpakDly;
  end
  n = n(:);
  % Make b and a rows
  b = b(:).';
  a = a(:).';
  
  % Make b and a of equal length
  [b,a]=eqtflength(b,a);
  
  % Compute group delay:
  [gd,w] = feval(dlyFcn, b,a,n,Fs,isNormalizedFreq,isWholeUnitCircle);
  
  % Compute frequency vector (normalized or  Hz):
  if isNormalizedFreq,
    f = w;
  else
    f = w * Fs/(2*pi);
  end  
else
  [gd, f] = grpdelay(b(1,1:3), b(1,4:6), n, inputArgs{:});
  for indx = 2:size(b, 1)
      gd = gd + grpdelay(b(indx,1:3), b(indx,4:6), n, inputArgs{:});    
  end    
end

if nargout == 0,
  % Produce plots of group delay calculations:
  newplot;
  if isNormalizedFreq,
    plot(f/pi,gd)
    xlabel(getString(message('signal:sigtools:siggui:NormalizedFrequencyRadsample','\times\pi')));
  else
    plot(f,gd)
    xlabel(getString(message('signal:grpdelay:FrequencyHz')))
  end
  ylabel(getString(message('signal:grpdelay:GroupDelaysamples')))
  set(gca,'xgrid','on','ygrid','on')
  
elseif nargout >= 1,  
  % Cast to enforce precision rules
  if isSingleArith
    gd_out = single(gd);
  else
    gd_out = gd;
  end
  if nargout == 2,
    % Cast to enforce precision rules
    if isSingleArith
      w_out = single(f);
    else
      w_out = f;
    end
  end
end

% ------------------------------------------------------------------
function [gd,w] = smithDly(b,a,n,Fs,isNormalizedFreq,isWholeUnitCircle)
%smithDly Computes group delay using the Smith algorithm

na = length(a);
c = conv(b, conj(a(na:-1:1)));
c = c(:).';	% make a row vector
nc = length(c);
cr = c.*(0:(nc-1));

if isWholeUnitCircle, s=1; else s=2; end

if length(n)==1,
   w = (2*pi/s*(0:n-1)/n)';
   if s*n >= nc	% pad with zeros to get the n values needed
      % dividenowarn temporarily suppresses warnings to avoid "Divide by zero"
      gd = dividenowarn(fft([cr zeros(1,s*n-nc)]),...
                        fft([c zeros(1,s*n-nc)]));
      gd = real(gd(1:n)) - ones(1,n)*(na-1);
   else	% find multiple of s*n points greater than nc
      nfact = s*ceil(nc/(s*n));
      mmax = n*nfact;
      % dividenowarn temporarily suppresses warnings to avoid "Divide by zero"
      gd = dividenowarn(fft(cr,mmax), fft(c,mmax));
      gd = real(gd(1:nfact:mmax)) - ones(1,n)*(na-1);
   end
   gd = gd(:);
else
    if isNormalizedFreq,
       w = n;
    else
       w = 2*pi*n/Fs;
    end
    s = exp(1i*w);
    gd = real(polyval(cr,s)./polyval(c,s));
    gd = gd - ones(size(gd))*(na-1);
end

% Linear phase FIRs
if signalpolyutils('islinphase',b,a,0),
    % Remove leading and trailing zeros
    startidx = find(b,1);
    stopidx = find(b,1,'last');
    if max(abs(b)) == 0,
        b = 0;
    else
        % Remove leading and trailing zeros of b
        b = b(startidx:stopidx);
    end

    % Compute group delay
    if isempty(startidx),
        G1 = 0;
    else
        G1 = max(startidx-1,0); % Delay introduced by leading zeros
    end
    G2 = (length(b)-1)/2; % Delay of symmetric FIR filter
    G = G1+G2;
    gd(:) = G;
end

gd = gd(:);

% ------------------------------------------------------------------
function [gd,w] = shpakDly(b,a,n,Fs,isNormalizedFreq,isWholeUnitCircle)
%shpakDly Computes group delay using the Shpak algorithm
% Tolerance for classifying roots as on the unit circle or at the origin
toler = 5*eps;

% Compute frequency vector (normalized or Hz):
if length(n) == 1
   if isWholeUnitCircle, ss=2*pi; else ss=pi; end
   w = (0:n-1)'/n * ss;
elseif isNormalizedFreq,
 	w = n;
else
 	w = 2*pi*n/Fs;
end

T = 1;
cw  = cos(w*T);
cw2 = cos(2*w*T);
sw  = sin(w*T);
sw2 = sin(2*w*T);
gd = zeros(length(w), 1);

if isreal(b) && isreal(a),
    % Get the real-valued second-order sections
    [sos,h0] = tf2sos(b,a); %#ok
    M = size(sos, 1);
else
    % Compute complex-valued second-order sections
    [z,p] = tf2zp(b,a);
    % Remove zeros and poles at the origin. 
    k = find(abs(z) > toler);
    % Each zero at the origin removes a unit delay
    gd = gd - T*(length(z)-length(k));
    z = z(k);
    k = find(abs(p) > toler);
    % Each pole at the origin adds a unit delay
    gd = gd + T*(length(p)-length(k));
    p = p(k);
    % Remove any zeros on the unit circle
    k = find(abs(abs(z) - 1) > toler);
    % Each zero on the unit circle adds a half-unit delay
    gd = gd + T*(length(z)-length(k)) / 2;
    z = z(k);
    % Again, but for poles
    k = find(abs(abs(p) - 1) > toler);
    gd = gd - T*(length(p)-length(k)) / 2;
    p = p(k);
    % Set up an array to hold the SOS
    numOrd = length(z);
    denOrd = length(p);
    M = max(ceil(numOrd/2), ceil(denOrd/2));
    sos = zeros(M,6);
    sos(:,1) = 1;
    sos(:,4) = 1;
    % Numerator
    sections = fix(numOrd/2);
    % No sorting, just pairing
    z1 = z(1:2:2*sections);
    z2 = z(2:2:2*sections);
    sos(1:sections,2:3) = [-(z1+z2) z1.*z2];
    if rem(numOrd, 2) == 1
        sos(sections+1,2) = -z(numOrd);
    end
    % Denominator
    sections = fix(denOrd/2);
    z1 = p(1:2:2*sections);
    z2 = p(2:2:2*sections);
    sos(1:sections,5:6) = [-(z1+z2) z1.*z2];
    if rem(denOrd, 2) == 1
        sos(sections+1,5) = -p(denOrd);
    end
end

for k=1:M
    % Numerator
    gd = gd + grpSection(sos(k,1:3),T,cw,sw,cw2,sw2);
    gd = gd - grpSection(sos(k,4:6),T,cw,sw,cw2,sw2);
end

% Compute the group delay for a first- or second-order polynomial
function gd=grpSection(q,T,cw,sw,cw2,sw2)
gd = zeros(length(cw), 1);
br = real(q);
bi = imag(q);
% Tolerance for finding symmetric and anti-symmetric polynomials
% (i.e., roots on the unit circle)
toler = 10*eps;

if q(2:3) == [0 0], %#ok
    % Zeroth-order section
    % (nothing to do!)
    
elseif q(3) == 0,
    % First-order section
    if isreal(q) && abs(br(1) - br(2)) < toler
        % Root at z=-1
        gd = T/2;
    elseif isreal(q) && abs(br(1) + br(2)) < toler
        % Root at z=1
        gd = T/2;
    else            
        b1 = br(1); b2 = br(2);
        g1 = bi(1); g2 = bi(2);
        u = g1*cw + b1*sw + g2;
        v = b1*cw - g1*sw + b2;
        du = T*(-g1*sw + b1*cw);
        dv = T*(-b1*sw - g1*cw);
    
        u2v2 = (b1^2 + g1^2 + b2^2 + g2^2) + 2*(b1*b2 + g1*g2)*cw + 2*(b1*g2 - b2*g1)*sw;
    
        % The following division could be zero/zero if we evaluate at a singularity
        k = find(abs(u2v2) > eps^(2/3));
        gd(k) = gd(k) + T - (v(k).*du(k) - u(k).*dv(k))./u2v2(k);
    end
else
    % Second-order section
    % First, check for symmetric and anti-symmetric polynomials
    if isreal(q) && all(abs(br - fliplr(br)) < toler)
        gd = T;
    elseif isreal(q) && all(abs(br + fliplr(br)) < toler)
        gd = T;
    else
        b1 = br(:,1); b2 = br(:,2); b3 = br(:,3);
        g1 = bi(:,1); g2 = bi(:,2); g3 = bi(:,3);
        u = g1*cw2 + b1*sw2 + g2*cw + b2*sw + g3;
        v = b1*cw2 - g1*sw2 + b2*cw - g2*sw + b3;
    
        du =  T*(-2*g1*sw2 + 2*b1*cw2 - g2*sw + b2*cw);
        dv = -T*( 2*b1*sw2 + 2*g1*cw2 + b2*sw + g2*cw);
    
        u2v2 = (b1^2 + g1^2) + (b2^2 + g2^2) + (b3^2 + g3^2) + ...
            2*(b1*b2 + g1*g2 + b2*b3 + g2*g3)*cw + 2*(b1*g2 - b2*g1 + b2*g3 - b3*g2)*sw + ...
            2*(b1*b3 + g1*g3)*cw2 + 2*(b1*g3 - b3*g1)*sw2; 
    
	    % The 2T compensates for using powers of +z rather than -z in the preceding derivatives. 
%         k = find(abs(u2v2) > eps^(2/3));
        k = 1:length(v);
        gd(k) = gd(k) + 2*T - (v(k).*du(k) - u(k).*dv(k))./u2v2(k);
    end
end   

% [EOF] grpdelay.m
