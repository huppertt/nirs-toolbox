function [h,a]=intfilt(R,L,freqmult)
%INTFILT Interpolation FIR Filter Design.
%   B = INTFILT(L,P,ALPHA) designs a linear phase FIR filter which, when 
%   used on a sequence interspersed with L-1 consecutive zeros every L
%   samples, performs bandlimited interpolation using the nearest 2*P
%   non-zero samples and assuming an original bandlimitedness of ALPHA
%   times the Nyquist frequency.  B is length 2*L*P-1.
%
%   The bandlimitedness factor expresses how much of the Nyquist interval
%   is occupied by the signal to be interpolated. This factor should be
%   greater than zero and less or equal to one. A bandlimitedness factor
%   that is less than one allows for a larger transition band and for
%   "don't-care" regions in the stopband of the filter where the signal has
%   no significant spectral content. The net result of specifying ALPHA
%   less than one is better stopband attenuation for a given L and P. See
%   example below.
%
%   B = INTFILT(L,N,'Lagrange') designs an FIR filter which performs Nth 
%   order Lagrange polynomial interpolation on a sequence that is inter-
%   spersed with L-1 consecutive zeros every L samples.  B is length 
%   (N+1)*L-1 for N odd and (N+1)*L for N even.  If both N and L are even, 
%   the filter designed is NOT linear phase.
%
%   Both types of filters are basically lowpass and have a gain of L in
%   the passband.
%
%   Example: Compare two filter with different bandlimitedness factors.
%   b1 = intfilt(5,10,1);  % Signal occupies full Nyquist interval
%   b2 = intfilt(5,10,.8); % Signal occupies 80% of Nyquist interval
%   hfvt = fvtool(b1,1,b2,1); legend(hfvt,'Alpha = 1','Alpha = 0.8');
%
%   See also INTERP, DECIMATE, RESAMPLE.

%   Reference:   Oetken, Parks, and Schuessler, "New Results in the
%      design of digital interpolators," IEEE Trans. Acoust., Speech,
%      Signal Processing, vol. ASSP-23, pp. 301-309, June 1975.

%   Copyright 1992-2013 The MathWorks, Inc.

narginchk(3,3)

% Cast to enforce Precision Rules
R = signal.internal.sigcasttofloat(R,'double','intfilt','L','allownumeric');
L = signal.internal.sigcasttofloat(L,'double','intfilt','P','allownumeric');                                 

if nargin == 3,
    if ischar(freqmult),
        type = freqmult;
        n = L;
    else
        % Cast to enforce Precision Rules
        freqmult = double(freqmult);
        type = 'b';
    end
end 

if strcmp(type(1),'b') || strcmp(type(1),'B'),
        
    n=2*R*L-1;

    if freqmult == 1
        
        M = [R R 0 0];
        F = [0 1/(2*R) 1/(2*R) .5];
    else
        
        M=R*[1 1];
        F=[0 freqmult/2/R];

        for f=(1/R):(1/R):.5,
            F=[F f-(freqmult/2/R) f+(freqmult/2/R)]; %#ok<*AGROW>
            M=[M 0 0];
        end;

        if (F(length(F))>.5),
            F(length(F))=.5;
        end;
    end

    h=firls(n-1,F*2,M);

elseif strcmp(type(1),'l') || strcmp(type(1),'L'),

% Inputs:
%    n   order of polynomial interpolation  (should be ODD!!!!)
%    R   discrete time sampling period (input to filter assumed 0 else)

    if n == 0
        h = ones(1,R);
        return
    end

    if ~(rem(n,2) || rem(R,2))
       warning(message('signal:intfilt:InvalidParamLinPhase', 'R', 'N'));
    end

    t=0:n*R+1;
    l=ones(n+1,length(t));
    for i=1:n+1
        for j=1:n+1
            if (j~=i)
                l(i,:)=l(i,:).*(t/R-j+1)/(i-j);
            end;
        end;
    end;

%    plot(t/R,l')
%    title('Lagrange Polynomials');
%    xlabel('time/R');
%    grid
%    l
%    pause;

    h=zeros(1,(n+1)*R);  

    for i=0:R-1
        for j=0:n
            h(j*R+i+1)=l((n-j)+1,round((n-1)/2*R+i+1));
        end;
    end;

    if h(1) == 0,
        h(1) = [];
    end

else
    error(message('signal:intfilt:InvalidParamFilterType'))
end

 if nargout > 1   
     a = 1;   
 end
