function [h,a]=firls(N,F,M,W,ftype)
% FIRLS Linear-phase FIR filter design using least-squares error minimization.
%   B=FIRLS(N,F,A) returns a length N+1 linear phase (real, symmetric
%   coefficients) FIR filter which has the best approximation to the
%   desired frequency response described by F and A in the least squares
%   sense. F is a vector of frequency band edges in pairs, in ascending
%   order between 0 and 1. 1 corresponds to the Nyquist frequency or half
%   the sampling frequency. A is a real vector the same size as F
%   which specifies the desired amplitude of the frequency response of the
%   resultant filter B. The desired response is the line connecting the
%   points (F(k),A(k)) and (F(k+1),A(k+1)) for odd k; FIRLS treats the
%   bands between F(k+1) and F(k+2) for odd k as "transition bands" or
%   "don't care" regions. Thus the desired amplitude is piecewise linear
%   with transition bands.  The integrated squared error is minimized.
%
%   For filters with a gain other than zero at Fs/2, e.g., highpass
%   and bandstop filters, N must be even.  Otherwise, N will be
%   incremented by one. Alternatively, you can use a trailing 'h' flag to
%   design a type 4 linear phase filter and avoid incrementing N.
%
%   B=FIRLS(N,F,A,W) uses the weights in W to weight the error. W has one
%   entry per band (so it is half the length of F and A) which tells
%   FIRLS how much emphasis to put on minimizing the integral squared error
%   in each band relative to the other bands.
%
%   B=FIRLS(N,F,A,'Hilbert') and B=FIRLS(N,F,A,W,'Hilbert') design filters
%   that have odd symmetry, that is, B(k) = -B(N+2-k) for k = 1, ..., N+1.
%   A special case is a Hilbert transformer which has an approx. amplitude
%   of 1 across the entire band, e.g. B=FIRLS(30,[.1 .9],[1 1],'Hilbert').
%
%   B=FIRLS(N,F,A,'differentiator') and B=FIRLS(N,F,A,W,'differentiator')
%   also design filters with odd symmetry, but with a special weighting
%   scheme for non-zero amplitude bands. The weight is assumed to be equal
%   to the inverse of frequency, squared, times the weight W. Thus the
%   filter has a much better fit at low frequency than at high frequency.
%   This designs FIR differentiators.
%
%   % Example of a length 31 lowpass filter.
%   h=firls(30,[0 .1 .2 .5]*2,[1 1 0 0]);
%   fvtool(h);
%
%   % Example of a length 45 lowpass differentiator.
%   h=firls(44,[0 .3 .4 1],[0 .2 0 0],'differentiator');
%   fvtool(h);
%
%   % Example of a length 26 type 4 highpass filter.
%   h=firls(25,[0 .4 .5 1],[0 0 1 1],'h');
%   fvtool(h);
%
%   See also FIRPM, FIR1, FIR2, FREQZ, FILTER, DESIGNFILT.

%       Author(s): T. Krauss
%   History: 10-18-91, original version
%            3-30-93, updated
%            9-1-95, optimize adjacent band case
%   Copyright 1988-2013 The MathWorks, Inc.
%     

% check number of arguments, set up defaults.
narginchk(3,5);

% Cast to enforce precision rules
N = signal.internal.sigcasttofloat(N,'double','firls','N','allownumeric');
F = signal.internal.sigcasttofloat(F,'double','firls','F','allownumeric');
M = signal.internal.sigcasttofloat(M,'double','firls','A','allownumeric');

if (max(F)>1) || (min(F)<0)
    error(message('signal:firls:InvalidRange'))
end
if (rem(length(F),2)~=0)
    error(message('signal:firls:MustHaveEvenLength', 'F'));
end
if (length(F) ~= length(M))
    error(message('signal:firls:UnequalLengths', 'F', 'A'));
end
if (nargin==3),
    W = ones(length(F)/2,1);
    ftype = '';
end
if (nargin==4),
    if ischar(W),
        ftype = W; 
        W = ones(length(F)/2,1);
    else
        ftype = '';
    end
end
if (nargin==5),
    if isempty(W),
        W = ones(length(F)/2,1);
    end
end
if isempty(ftype)
    ftype = 0;  differ = 0;
else
    ftype = lower(ftype);
    if strcmpi(ftype,'h') || strcmpi(ftype,'hilbert')
        ftype = 1;  differ = 0;
    elseif strcmpi(ftype,'d') || strcmpi(ftype,'differentiator')
        ftype = 1;  differ = 1;
    else
        error(message('signal:firls:InvalidEnum'))
    end
end
% Cast to enforce precision rules
W = signal.internal.sigcasttofloat(W,'double','firls','W','allownumeric');

% Check for valid filter length
[N,msg1,msg2,msgobj] = firchk(N,F(end),M,ftype);
if ~isempty(msg1), error(msgobj); end;
if ~isempty(msg2), warning(msgobj); end;

N = N+1;                   % filter length
F=F(:)/2;  M=M(:);  W=sqrt(W(:));  % make these guys columns
dF = diff(F);

if (length(F) ~= length(W)*2)
    error(message('signal:firls:InvalidDimensions'));
end
if any(dF<0),
    error(message('signal:firls:InvalidFreqVec'))
end

% Fix for 67187
if all(dF(2:2:length(dF)-1)==0) && length(dF) > 1,
    fullband = 1;
else
    fullband = 0;
end
if all((W-W(1))==0)
    constant_weights = 1;
else
    constant_weights = 0;
end

L=(N-1)/2;

Nodd = rem(N,2);

if (ftype == 0),  % Type I and Type II linear phase FIR
    % basis vectors are cos(2*pi*m*f) (see m below)
    if ~Nodd
        m=(0:L)+.5;   % type II
    else
        m=(0:L);      % type I
    end
    k=m';
    need_matrix = (~fullband) || (~constant_weights);
    if need_matrix
        I1=k(:,ones(size(m)))+m(ones(size(k)),:);    % entries are m + k
        I2=k(:,ones(size(m)))-m(ones(size(k)),:);    % entries are m - k
        G=zeros(size(I1));
    end

    if Nodd
        k=k(2:length(k));
        b0=0;       %  first entry must be handled separately (where k(1)=0)
    end;
    b=zeros(size(k));
    for s=1:2:length(F),
        m=(M(s+1)-M(s))/(F(s+1)-F(s));    %  slope
        b1=M(s)-m*F(s);                   %  y-intercept
        if Nodd
            b0 = b0 + (b1*(F(s+1)-F(s)) + m/2*(F(s+1)*F(s+1)-F(s)*F(s)))...
                * abs(W((s+1)/2)^2) ;
        end
        b = b+(m/(4*pi*pi)*(cos(2*pi*k*F(s+1))-cos(2*pi*k*F(s)))./(k.*k))...
            * abs(W((s+1)/2)^2);
        b = b + (F(s+1)*(m*F(s+1)+b1)*sinc(2*k*F(s+1)) ...
            - F(s)*(m*F(s)+b1)*sinc(2*k*F(s))) ...
            * abs(W((s+1)/2)^2);
        if need_matrix
            G = G + (.5*F(s+1)*(sinc(2*I1*F(s+1))+sinc(2*I2*F(s+1))) ...
                - .5*F(s)*(sinc(2*I1*F(s))+sinc(2*I2*F(s))) ) ...
                * abs(W((s+1)/2)^2);
        end
    end;
    if Nodd
        b=[b0; b];
    end;

    if need_matrix
        a=G\b;
    else
        a=(W(1)^2)*4*b;
        if Nodd
            a(1) = a(1)/2;
        end
    end
    if Nodd
        h=[a(L+1:-1:2)/2; a(1); a(2:L+1)/2].';
    else
        h=.5*[flipud(a); a].';
    end;
elseif (ftype == 1),  % Type III and Type IV linear phase FIR
    %  basis vectors are sin(2*pi*m*f) (see m below)
    if (differ),      % weight non-zero bands with 1/f^2
        do_weight = ( abs(M(1:2:length(M))) +  abs(M(2:2:length(M))) ) > 0;
    else
        do_weight = zeros(size(F));
    end

    if Nodd
        m=(1:L);      % type III
    else
        m=(0:L)+.5;   % type IV
    end;
    k=m';
    b=zeros(size(k));

    need_matrix = (~fullband) || (any(do_weight)) || (~constant_weights);
    if need_matrix
        I1=k(:,ones(size(m)))+m(ones(size(k)),:);    % entries are m + k
        I2=k(:,ones(size(m)))-m(ones(size(k)),:);    % entries are m - k
        G=zeros(size(I1));
    end

    for s=1:2:length(F),
        if (do_weight((s+1)/2)),      % weight bands with 1/f^2
            if F(s) == 0, F(s) = 1e-5; end     % avoid singularities
            m=(M(s+1)-M(s))/(F(s+1)-F(s));
            b1=M(s)-m*F(s);
            snint1 = sineint(2*pi*k*F(s+1)) - sineint(2*pi*k*F(s));
            %snint1 = (-1/2/i)*(expint(1i*2*pi*k*F(s+1)) ...
            %    -expint(-i*2*pi*k*F(s+1)) -expint(1i*2*pi*k*F(s)) ...
            %    +expint(-i*2*pi*k*F(s)) );
            % csint1 = cosint(2*pi*k*F(s+1)) - cosint(2*pi*k*F(s)) ;
            csint1 = (-1/2)*(expint(1i*2*pi*k*F(s+1))+expint(-1i*2*pi*k*F(s+1))...
                -expint(1i*2*pi*k*F(s))  -expint(-1i*2*pi*k*F(s)) );
            b=b + ( m*snint1 ...
                + b1*2*pi*k.*( -sinc(2*k*F(s+1)) + sinc(2*k*F(s)) + csint1 ))...
                * abs(W((s+1)/2)^2);
            snint1 = sineint(2*pi*F(s+1)*(-I2));
            snint2 = sineint(2*pi*F(s+1)*I1);
            snint3 = sineint(2*pi*F(s)*(-I2));
            snint4 = sineint(2*pi*F(s)*I1);
            G = G - ( ( -1/2*( cos(2*pi*F(s+1)*(-I2))/F(s+1)  ...
                - 2*snint1*pi.*I2 ...
                - cos(2*pi*F(s+1)*I1)/F(s+1) ...
                - 2*snint2*pi.*I1 )) ...
                - ( -1/2*( cos(2*pi*F(s)*(-I2))/F(s)  ...
                - 2*snint3*pi.*I2 ...
                - cos(2*pi*F(s)*I1)/F(s) ...
                - 2*snint4*pi.*I1) ) ) ...
                * abs(W((s+1)/2)^2);
        else      % use usual weights
            m=(M(s+1)-M(s))/(F(s+1)-F(s));
            b1=M(s)-m*F(s);
            b=b+(m/(4*pi*pi)*(sin(2*pi*k*F(s+1))-sin(2*pi*k*F(s)))./(k.*k))...
                * abs(W((s+1)/2)^2) ;
            b = b + (((m*F(s)+b1)*cos(2*pi*k*F(s)) - ...
                (m*F(s+1)+b1)*cos(2*pi*k*F(s+1)))./(2*pi*k)) ...
                * abs(W((s+1)/2)^2) ;
            if need_matrix
                G = G + (.5*F(s+1)*(sinc(2*I1*F(s+1))-sinc(2*I2*F(s+1))) ...
                    - .5*F(s)*(sinc(2*I1*F(s))-sinc(2*I2*F(s)))) * ...
                    abs(W((s+1)/2)^2);
            end
        end;
    end

    if need_matrix
        a=G\b;
    else
        a=-4*b*(W(1)^2);
    end
    if Nodd
        h=.5*[flipud(a); 0; -a].';
    else
        h=.5*[flipud(a); -a].';
    end
    if differ, h=-h; end
end

if nargout > 1
    a = 1;
end

%----------------------------------------------------------------------------
function y = sineint(x)
% SINEINT (a.k.a. SININT)   Numerical Sine Integral
%   Used by FIRLS in the Signal Processing Toolbox.
%   Untested for complex or imaginary inputs.
%
%   See also SININT in the Symbolic Toolbox.

%   Was Revision: 1.5, Date: 1996/03/15 20:55:51

i1 = find(real(x)<0);   % this equation is not valid if x is in the
% left-hand plane of the complex plane.
% use relation Si(-z) = -Si(z) in this case (Eq 5.2.19, Abramowitz
%  & Stegun).
x(i1) = -x(i1);
y = zeros(size(x));
ind = find(x);
% equation 5.2.21 Abramowitz & Stegun
%  y(ind) = (1/(2*i))*(expint(1i*x(ind)) - expint(-i*x(ind))) + pi/2;
y(ind) = imag(expint(1i*x(ind))) + pi/2;
y(i1) = -y(i1);

