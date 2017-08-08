function [zd, pd, kd, dd] = bilinear(z, p, k, fs, fp, fp1)
%BILINEAR Bilinear transformation with optional frequency prewarping.
%   [Zd,Pd,Kd] = BILINEAR(Z,P,K,Fs) converts the s-domain transfer
%   function specified by Z, P, and K to a z-transform discrete
%   equivalent obtained from the bilinear transformation:
%
%      H(z) = H(s) |
%                  | s = 2*Fs*(z-1)/(z+1)
%
%   where column vectors Z and P specify the zeros and poles, scalar
%   K specifies the gain, and Fs is the sample frequency in Hz.
%
%   [NUMd,DENd] = BILINEAR(NUM,DEN,Fs), where NUM and DEN are
%   row vectors containing numerator and denominator transfer
%   function coefficients, NUM(s)/DEN(s), in descending powers of
%   s, transforms to z-transform coefficients NUMd(z)/DENd(z).
%
%   [Ad,Bd,Cd,Dd] = BILINEAR(A,B,C,D,Fs) is a state-space version.
%
%   Each of the above three forms of BILINEAR accepts an optional
%   additional input argument that specifies prewarping.
%
%   For example, [Zd,Pd,Kd] = BILINEAR(Z,P,K,Fs,Fp) applies prewarping
%   before the bilinear transformation so that the frequency responses
%   before and after mapping match exactly at frequency point Fp
%   (match point Fp is specified in Hz).
%
%   % Example:
%   %   Design a 6-th order Elliptic analog low pass filter and transform
%   %   it to a Discrete-time representation.
%
%   Fs =0.5;                            % Sampling Frequency
%   [z,p,k]=ellipap(6,5,90);            % Lowpass filter prototype
%   [num,den]=zp2tf(z,p,k);             % Convert to transfer function form
%   [numd,dend]=bilinear(num,den,Fs);   % Analog to Digital conversion
%   fvtool(numd,dend)                   % Visualize the filter
%
%   See also IMPINVAR.

%   Author(s): J.N. Little, 4-28-87
%   	   J.N. Little, 5-5-87, revised
%   Copyright 1988-2006 The MathWorks, Inc.

%   Gene Franklin, Stanford Univ., motivated the state-space
%   approach to the bilinear transformation.

[mn,nn] = size(z);
[md,nd] = size(p);
if (nd == 1 && nn < 2) && nargout ~= 4	% In zero-pole-gain form
    if mn > md
        error(message('signal:bilinear:InvalidRange'))
    end
    if nargin == 5
      % Cast to enforce Precision Rules
      fp = signal.internal.sigcasttofloat(fp,'double','bilinear','',...
        'allownumeric');
      fs = signal.internal.sigcasttofloat(fs,'double','bilinear','',...
        'allownumeric');
      % Prewarp
      fp = 2*pi*fp;
      fs = fp/tan(fp/fs/2);
    else
      % Cast to enforce Precision Rules
      fs = 2*signal.internal.sigcasttofloat(fs,'double','bilinear','',...
        'allownumeric');
    end
    % Cast to enforce Precision Rules
    if any([signal.internal.sigcheckfloattype(z,'single','bilinear','Z') ...
        signal.internal.sigcheckfloattype(p,'single','bilinear','P')...
        signal.internal.sigcheckfloattype(k,'single','bilinear','K')])
      z = single(z);
      p = single(p);
      k = single(k);
    end
    z = z(isfinite(z));	 % Strip infinities from zeros
    pd = (1+p/fs)./(1-p/fs); % Do bilinear transformation
    zd = (1+z/fs)./(1-z/fs);
    % real(kd) or just kd?
    kd = (k*prod(fs-z)./prod(fs-p));
    zd = [zd;-ones(length(pd)-length(zd),1)];  % Add extra zeros at -1
    
elseif (md == 1 && mn == 1) || nargout == 4 %
    if nargout == 4		% State-space case
        % Cast to enforce Precision Rules
        a = z;
        b = p;
        c = k;
        d = fs; 
        fs = signal.internal.sigcasttofloat(fp,'double','bilinear','',...
          'allownumeric');
        if any([signal.internal.sigcheckfloattype(a,'single','bilinear','A')...
            signal.internal.sigcheckfloattype(b,'single','bilinear','B')...
            signal.internal.sigcheckfloattype(c,'single','bilinear','C')...
            signal.internal.sigcheckfloattype(d,'single','bilinear','D')])
          a = single(a);
          b = single(b);
          c = single(c);
          d = single(d);
        end
        error(abcdchk(a,b,c,d));
        if nargin == 6			% Prewarp
          % Cast to enforce Precision Rules
          fp = signal.internal.sigcasttofloat(fp1,'double','bilinear','',...
            'allownumeric');		% Decode arguments
          fp = 2*pi*fp;
          fs = fp/tan(fp/fs/2)/2;
        end
    else			% Transfer function case
        if nn > nd
            error(message('signal:bilinear:InvalidRange'))
        end
        num = z; den = p;		% Decode arguments
        if any([signal.internal.sigcheckfloattype(num,'single','bilinear','NUM')...
            signal.internal.sigcheckfloattype(den,'single','bilinear','DEN')])
          num = single(num);
          den = single(den);
        end
        if nargin == 4			% Prewarp
            % Cast to enforce Precision Rules
            fp = signal.internal.sigcasttofloat(fs,'double','bilinear','',...
              'allownumeric'); 
            fs = signal.internal.sigcasttofloat(k,'double','bilinear','',...
              'allownumeric');	% Decode arguments
            fp = 2*pi*fp;
            fs = fp/tan(fp/fs/2)/2;
            
        else
            % Cast to enforce Precision Rules
            fs = signal.internal.sigcasttofloat(k,'double','bilinear','',...
              'allownumeric');	% Decode arguments
        end
        
        % Put num(s)/den(s) in state-space canonical form.
        [a,b,c,d] = tf2ss(num,den);
    end
    % Now do state-space version of bilinear transformation:
    t = 1/fs;
    r = sqrt(t);
    t1 = eye(size(a)) + a*t/2;
    t2 = eye(size(a)) - a*t/2;
    ad = t2\t1;
    bd = t/r*(t2\b);
    cd = r*c/t2;
    dd = c/t2*b*t/2 + d;
    if nargout == 4
        zd = ad; pd = bd; kd = cd;
    else
        % Convert back to transfer function form:
        p = poly(ad);
        zd = poly(ad-bd*cd)+(dd-1)*p;
        pd = p;
    end
else
    error(message('signal:bilinear:SignalErr'))
end
